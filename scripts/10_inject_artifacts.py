import argparse, gzip, json, os, shutil, sys
import numpy as np, pandas as pd
import pyarrow.parquet as pq
import scipy.io, scipy.sparse as sp

BIN_SIZES = [8, 16, 32, 64]
BASE_UM = 8

def log(*a): print(*a, flush=True)

def read_mex(d):
    m = scipy.io.mmread(os.path.join(d, "matrix.mtx.gz")).tocsc()
    with gzip.open(os.path.join(d, "barcodes.tsv.gz"), "rt") as f:
        bc = np.array([l.strip() for l in f])
    feat = pd.read_csv(os.path.join(d, "features.tsv.gz"), sep="\t", header=None,
                       names=["id", "symbol", "type"])
    return m, bc, feat

def write_mex(out, mat, barcodes, feat):
    os.makedirs(out, exist_ok=True)
    with gzip.open(os.path.join(out, "matrix.mtx.gz"), "wb") as f:
        scipy.io.mmwrite(f, mat.tocoo(), field="integer")
    with gzip.open(os.path.join(out, "barcodes.tsv.gz"), "wt") as f:
        f.write("\n".join(barcodes) + "\n")
    feat.to_csv(os.path.join(out, "features.tsv.gz"), sep="\t",
                header=False, index=False, compression="gzip")

def grid_basis(pos):
    A = np.column_stack([pos.array_row.to_numpy(float),
                         pos.array_col.to_numpy(float), np.ones(len(pos))])
    B = pos[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy(float)
    coef, *_ = np.linalg.lstsq(A, B, rcond=None)
    resid = np.abs(B - A @ coef).max()
    return coef, resid

def truncated_mask(meta, arr_max_row, arr_max_col, f):
    return (meta.array_row.to_numpy() == arr_max_row // f) | \
           (meta.array_col.to_numpy() == arr_max_col // f)

def pick_regions(pos8, diam_um, n_regions, rng, edge_margin_bins, arr_max,
                 location="interior"):
    r_bins = (diam_um / 2.0) / BASE_UM
    tis = pos8[pos8.in_tissue == 1]
    rmax, cmax = arr_max
    m = edge_margin_bins + r_bins
    cand = tis[(tis.array_row > m) & (tis.array_row < rmax - m) &
               (tis.array_col > m) & (tis.array_col < cmax - m)]
    if location == "edge":
        from scipy.ndimage import distance_transform_edt
        grid = np.zeros((int(rmax) + 1, int(cmax) + 1), bool)
        grid[pos8.array_row.to_numpy()[pos8.in_tissue.to_numpy() == 1],
             pos8.array_col.to_numpy()[pos8.in_tissue.to_numpy() == 1]] = True
        dist = distance_transform_edt(grid)
        d_at = dist[cand.array_row.to_numpy(), cand.array_col.to_numpy()]
        cand = cand[(d_at > 0) & (d_at <= r_bins)]
    if len(cand) == 0:
        return []
    order = rng.permutation(len(cand))
    centres, min_sep = [], 2.0 * r_bins + 2.0
    for i in order:
        r, c = int(cand.array_row.iloc[i]), int(cand.array_col.iloc[i])
        if all((r - pr) ** 2 + (c - pc) ** 2 > min_sep ** 2 for pr, pc in centres):
            centres.append((r, c))
        if len(centres) >= n_regions:
            break
    return centres

def region_mask(pos8, centres, diam_um):
    r_bins = (diam_um / 2.0) / BASE_UM
    ar = pos8.array_row.to_numpy(); ac = pos8.array_col.to_numpy()
    mask = np.zeros(len(pos8), bool)
    for pr, pc in centres:
        mask |= ((ar - pr) ** 2 + (ac - pc) ** 2) <= r_bins ** 2
    return mask & (pos8.in_tissue.to_numpy() == 1)

def inject_dropout(mat, col_idx, keep_p, rng):
    out = mat.copy().tocsc()
    for j in col_idx:
        s, e = out.indptr[j], out.indptr[j + 1]
        if e > s:
            out.data[s:e] = rng.binomial(out.data[s:e].astype(np.int64), keep_p)
    out.eliminate_zeros()
    return out

def inject_mito(mat, col_idx, mito_rows, target_frac, rng):
    out = mat.copy().tocsc().astype(np.int64)
    mito_rows = np.asarray(mito_rows)
    weights = np.asarray(mat[mito_rows].sum(1)).ravel().astype(float)
    weights = (weights / weights.sum() if weights.sum() > 0
               else np.full(len(mito_rows), 1.0 / len(mito_rows)))

    add_rows, add_cols, add_vals = [], [], []
    for j in col_idx:
        s_, e_ = out.indptr[j], out.indptr[j + 1]
        total = out.data[s_:e_].sum()
        cur_mito = out.data[s_:e_][np.isin(out.indices[s_:e_], mito_rows)].sum()
        other = total - cur_mito
        need = target_frac * other / (1.0 - target_frac) - cur_mito
        if need <= 0:
            continue
        n_add = rng.poisson(need)
        if n_add == 0:
            continue
        counts = rng.multinomial(n_add, weights)
        nz = np.where(counts > 0)[0]
        add_rows.extend(mito_rows[nz]); add_cols.extend([j] * len(nz))
        add_vals.extend(counts[nz])

    if add_vals:
        add = sp.csc_matrix((add_vals, (add_rows, add_cols)), shape=out.shape)
        out = (out + add).tocsc()
    out.eliminate_zeros()
    return out

def aggregate(mat, pos8, f):
    pr = pos8.array_row.to_numpy() // f
    pc = pos8.array_col.to_numpy() // f
    key = pr.astype(np.int64) * (1 << 32) + pc
    uniq, inv = np.unique(key, return_inverse=True)
    M = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)),
                      shape=(len(inv), len(uniq)))
    agg = (mat.tocsr() @ M).tocsc()
    out_um = BASE_UM * f
    bc = [f"s_{out_um:03d}um_{int(k >> 32):05d}_{int(k & 0xFFFFFFFF):05d}-1" for k in uniq]
    meta = pd.DataFrame({"barcode": bc,
                         "array_row": (uniq >> 32).astype(int),
                         "array_col": (uniq & 0xFFFFFFFF).astype(int),
                         "n_children": np.bincount(inv)})
    return agg, meta

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, help="the section's binned_outputs directory")
    ap.add_argument("--out", required=True, help="the run directory to write the conditions into")
    ap.add_argument("--diameters", type=int, nargs="+", default=[16, 32, 64, 128, 256])
    ap.add_argument("--n-regions", type=int, default=25)
    ap.add_argument("--keep-p", type=float, default=0.2, help="dropout: fraction of counts retained")
    ap.add_argument("--mito-target-frac", type=float, default=0.25,
                help="mito: mitochondrial fraction to raise the region to")
    ap.add_argument("--mito-odds-ratio", type=float, default=None,
                help="mito: instead of an absolute target, raise the mitochondrial "
                     "odds of the region to this multiple of the section's own "
                     "median. Necrosis raises the mitochondrial fraction relative "
                     "to where the tissue started, so an absolute target is a "
                     "different insult in different tissue - in human kidney, whose "
                     "healthy bins sit at 20.7%, a 0.25 target is barely a "
                     "perturbation at all.")
    ap.add_argument("--mito-baseline", default="median", choices=("median", "native-mean"),
                    help="the section baseline --mito-odds-ratio multiplies. 'median' is the "
                         "median per-bin mitochondrial fraction of the thinned run, which is "
                         "zero at low depth because most bins then hold no mitochondrial "
                         "count, so the target moved with depth. 'native-mean' is the mean "
                         "per-bin fraction of the unthinned section, which thinning leaves "
                         "unchanged in expectation, so all depths of a section share one "
                         "target")
    ap.add_argument("--record-suffix", default="",
                    help="suffix for manifest.csv and injection_params.json, so that an "
                         "injection into an existing run keeps the earlier record")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--subsample", type=float, default=1.0,
                    help="binomially thin every count before injecting, to emulate a "
                         "shallower library. 1.0 leaves the matrix untouched.")
    ap.add_argument("--target-umi-per-bin", type=float, default=None,
                    help="thin to this mean UMI per bin instead of a fixed fraction. "
                         "Lets several datasets of different native depth be brought "
                         "to a common depth, which is what makes depth a controlled "
                         "variable rather than a property of the tissue.")
    ap.add_argument("--kinds", nargs="+", default=["dropout", "mito"],
                    help="which injection types to generate")
    ap.add_argument("--with-clean", action="store_true", default=False,
                    help="also write the uncorrupted condition")
    ap.add_argument("--location", default="interior", choices=("interior", "edge"),
                    help="where the lesions sit relative to the tissue boundary; "
                         "edge appends 'edge' to the condition label, with no underscore, "
                         "because the scorer reads the diameter from the first underscore")
    ap.add_argument("--kind-label", default=None,
                    help="name used in the condition tag; lets a second severity "
                         "live alongside the first (e.g. dropout05_64um)")
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    d8 = f"{a.src}/square_008um"
    mat, bc8, feat = read_mex(f"{d8}/filtered_feature_bc_matrix")
    pos8 = pq.read_table(f"{d8}/spatial/tissue_positions.parquet").to_pandas()
    sf8 = json.load(open(f"{d8}/spatial/scalefactors_json.json"))

    if a.mito_baseline == "native-mean":
        mr = np.where(feat.symbol.astype(str).str.upper().str.startswith("MT-"))[0]
        m0 = np.asarray(mat[mr].sum(0)).ravel(); t0 = np.asarray(mat.sum(0)).ravel()
        a.mito_baseline_value = float(np.mean(m0[t0 > 0] / t0[t0 > 0]))
        log(f"native mean per-bin mitochondrial fraction {a.mito_baseline_value:.4f}")
    if a.target_umi_per_bin is not None:
        native = mat.sum() / mat.shape[1]
        if a.target_umi_per_bin > native:
            sys.exit(f"target {a.target_umi_per_bin} UMIs/bin exceeds this sample's "
                     f"native depth {native:.1f}; counts can be thinned, not invented")
        a.subsample = a.target_umi_per_bin / native
        log(f"native {native:.1f} UMIs/bin -> target {a.target_umi_per_bin:.1f} "
            f"(fraction {a.subsample:.4f})")
    if a.subsample < 1.0:
        before = mat.sum() / mat.shape[1]
        mat = mat.tocsc()
        mat.data = rng.binomial(mat.data.astype(np.int64), a.subsample).astype(mat.dtype)
        mat.eliminate_zeros()
        log(f"subsampled to {a.subsample:.0%}: {before:.1f} -> "
            f"{mat.sum() / mat.shape[1]:.1f} UMIs/bin")

    pos8 = pos8.set_index("barcode").reindex(bc8).reset_index()
    assert pos8.barcode.tolist() == list(bc8)
    log(f"8 um matrix: {mat.shape[0]:,} genes x {mat.shape[1]:,} bins, "
        f"{mat.nnz:,} nonzero, mean {mat.sum() / mat.shape[1]:.1f} UMIs/bin")

    mito_rows = np.where(feat.symbol.astype(str).str.upper().str.startswith("MT-"))[0]
    log(f"mitochondrial probes: {len(mito_rows)}  ({', '.join(feat.symbol.iloc[mito_rows][:6])} ...)")
    base_mito = np.asarray(mat[mito_rows].sum(0)).ravel()
    tot = np.asarray(mat.sum(0)).ravel()
    log(f"baseline mito ratio: mean {np.mean(base_mito / np.maximum(tot, 1)):.4f}, "
        f"bins with zero mito counts: {np.mean(base_mito == 0) * 100:.1f}%")
    if a.mito_odds_ratio is not None:
        ok = tot > 0
        f0 = (a.mito_baseline_value if a.mito_baseline == "native-mean"
              else float(np.median(base_mito[ok] / tot[ok])))
        odds = a.mito_odds_ratio * f0 / (1.0 - f0)
        rel = odds / (1.0 + odds)
        a.mito_target_frac = max(a.mito_target_frac, rel)
        log(f"mito odds ratio {a.mito_odds_ratio:g} x section {a.mito_baseline} {f0:.4f} "
            f"-> relative target {rel:.4f}; using max(absolute, relative) = "
            f"{a.mito_target_frac:.4f}")

    full_pos8 = pq.read_table(f"{d8}/spatial/tissue_positions.parquet").to_pandas()
    coef, resid = grid_basis(full_pos8)
    log(f"grid basis residual: {resid:.4f} px")

    os.makedirs(a.out, exist_ok=True)
    manifest = []

    arr_max = (int(full_pos8.array_row.max()), int(full_pos8.array_col.max()))

    regions = {d: pick_regions(pos8, d, a.n_regions, rng, 8, arr_max, a.location)
               for d in a.diameters}
    for d, cs in regions.items():
        log(f"diameter {d:3d} um: {len(cs)} {a.location} regions")

    conditions = [("clean", None)] if a.with_clean else []
    for d in a.diameters:
        conditions += [(k, d) for k in a.kinds]

    seen_tags = set()
    for kind, diam in conditions:
        label = kind
        if a.kind_label:
            if "=" in a.kind_label:
                k_, l_ = a.kind_label.split("=", 1)
                if k_ == kind:
                    label = l_
            elif kind == "dropout":
                label = a.kind_label
        if a.location == "edge" and kind != "clean":
            label = f"{label}edge"
        tag = "clean" if kind == "clean" else f"{label}_{diam}um"
        if tag in seen_tags:
            sys.exit(f"condition tag collision: {tag} would be written twice")
        seen_tags.add(tag)
        log(f"\n=== {tag} ===")
        if kind == "clean":
            truth = np.zeros(len(bc8), bool); cur = mat; centres = []
        else:
            centres = regions[diam]
            truth = region_mask(pos8, centres, diam)
            cols = np.where(truth)[0]
            log(f"  {len(centres)} regions, {truth.sum():,} corrupted bins "
                f"({truth.mean() * 100:.2f}% of in-tissue)")
            cur = (inject_dropout(mat, cols, a.keep_p, rng) if kind == "dropout"
                   else inject_mito(mat, cols, mito_rows, a.mito_target_frac, rng))
            nm = np.asarray(cur[mito_rows].sum(0)).ravel()
            nt = np.asarray(cur.sum(0)).ravel()
            log(f"  library size in region : {tot[cols].mean():7.1f} -> {nt[cols].mean():7.1f}")
            log(f"  mito ratio in region   : {np.mean(base_mito[cols]/np.maximum(tot[cols],1)):.4f}"
                f" -> {np.mean(nm[cols]/np.maximum(nt[cols],1)):.4f}")

        root = os.path.join(a.out, tag)
        os.makedirs(root, exist_ok=True)
        pd.DataFrame({"barcode": bc8, "is_artifact": truth}).to_parquet(
            os.path.join(root, "ground_truth_008um.parquet"), index=False)

        for b in BIN_SIZES:
            f = b // BASE_UM
            outd = os.path.join(root, f"square_{b:03d}um")
            if f == 1:
                m_out, meta = cur, pos8[["barcode", "array_row", "array_col"]].copy()
                meta["n_children"] = 1
            else:
                m_out, meta = aggregate(cur, pos8, f)
            A = np.column_stack([meta.array_row.to_numpy(float),
                                 meta.array_col.to_numpy(float), np.ones(len(meta))])
            xy = A @ coef
            meta["pxl_row_in_fullres"] = xy[:, 0]
            meta["pxl_col_in_fullres"] = xy[:, 1]
            meta["in_tissue"] = 1
            meta["is_truncated"] = truncated_mask(meta, arr_max[0], arr_max[1], f)
            write_mex(os.path.join(outd, "filtered_feature_bc_matrix"), m_out, meta.barcode, feat)
            os.makedirs(os.path.join(outd, "spatial"), exist_ok=True)
            meta.to_parquet(os.path.join(outd, "spatial", "tissue_positions.parquet"), index=False)
            sf = dict(sf8); sf["bin_size_um"] = float(b)
            sf["spot_diameter_fullres"] = sf8["spot_diameter_fullres"] * f
            json.dump(sf, open(os.path.join(outd, "spatial", "scalefactors_json.json"), "w"), indent=4)
            manifest.append({"condition": tag, "kind": kind, "diameter_um": diam,
                             "bin_um": b, "n_bins": int(m_out.shape[1]),
                             "n_truncated": int(meta.is_truncated.sum()),
                             "path": outd})
            log(f"  {b:3d} um: {m_out.shape[1]:>7,} bins  "
                f"({int(meta.is_truncated.sum())} truncated)")

    sfx = a.record_suffix
    pd.DataFrame(manifest).to_csv(os.path.join(a.out, f"manifest{sfx}.csv"), index=False)
    json.dump(vars(a), open(os.path.join(a.out, f"injection_params{sfx}.json"), "w"), indent=2)
    log(f"\nwrote {a.out}/manifest{sfx}.csv")

if __name__ == "__main__":
    main()
