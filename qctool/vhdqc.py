import argparse, gzip, importlib.util, json, os, shutil, subprocess, sys, time
import numpy as np, pandas as pd, pyarrow.parquet as pq

HELP = """QC of one Visium HD section by the procedure of Fig. 6.

    python3 qctool/vhdqc.py --input .../outs/binned_outputs --out qc_out --damage-um 200
    python3 qctool/vhdqc.py --input ... --out qc_out --damage-um 60 --edge
    python3 qctool/vhdqc.py --input ... --out qc_out

Writes removed_8um.txt, bins_8um.csv.gz, qc_report.json, qc_report.md and qc_map.png to --out."""

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("inj", os.path.join(ROOT, "scripts", "10_inject_artifacts.py"))
inj = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inj)
EDGE_RADIUS_BINS = 16

def settings_path(a):
    for p in (a.settings, os.path.join(HERE, "guideline_settings.json"),
              os.path.join(ROOT, "benchmark_grid", "figures", "guideline_settings.json")):
        if p and os.path.exists(p):
            return p
    sys.exit("no guideline_settings.json: run scripts/128_guideline_figure.py or pass --settings")

def check_core():
    core = os.path.join(HERE, "qc_core.R")
    srcs = [os.path.join(ROOT, "scripts", "12_run_arms_k48all.R"), os.path.join(ROOT, "R", "ad.R")]
    if not os.path.exists(core) or any(os.path.exists(s) and os.path.getmtime(s) > os.path.getmtime(core)
                                       for s in srcs):
        sys.exit("qctool/qc_core.R is missing or older than its sources: run Rscript qctool/make_core.R")

def plan(S, mean_umi, damage_um, edge):
    tier = "hi" if mean_umi >= S["depth_cut_umi_per_bin"] else "lo"
    cut = S["depth_cut_umi_per_bin"]
    depth = f"{mean_umi:.0f} UMIs per 8 um bin ({'at or above' if tier == 'hi' else 'below'} {cut})"
    tests, notes = [], []
    if damage_um is None:
        b = S["mad"][tier]
        tests.append(dict(test="mad", bin=b, k=None,
                          why=f"no damage size given: section-wide MAD threshold on {b} um bins; {depth}"))
        if edge:
            notes.append("--edge is used only when a damage size is given (Fig. 6); ignored here")
    else:
        cls = next(c for c in S["classes"] if c["min_um"] <= damage_um and
                   (c["max_um"] is None or damage_um < c["max_um"]))
        r = cls[tier]
        tests.append(dict(test="spotsweeper", bin=r["bin"], k=r["k"],
                          why=f"damage of {damage_um:g} um is in the {cls['label']} class; {depth}"))
        if r.get("prov"):
            notes.append("this setting is provisional in Fig. 6 (fewer sections than required)")
        if edge:
            b = S["spatialartifacts"][tier]
            tests.append(dict(test="spatialartifacts", bin=b, k=None,
                              why=f"damage at the tissue edge: SpatialArtifacts on {b} um bins; {depth}"))
    return tier, tests, notes

def prepare(src, views, bins, mat, bc8, feat):
    d8 = os.path.join(src, "square_008um")
    full = pq.read_table(os.path.join(d8, "spatial", "tissue_positions.parquet")).to_pandas()
    pos8 = full.set_index("barcode").reindex(bc8).reset_index()
    assert pos8.barcode.tolist() == list(bc8)
    if pos8.array_row.isna().any():
        sys.exit("some filtered barcodes have no row in tissue_positions.parquet")
    sf8 = json.load(open(os.path.join(d8, "spatial", "scalefactors_json.json")))
    coef, resid = inj.grid_basis(full)
    if resid > 1.0:
        sys.exit(f"grid basis residual {resid:.3f} px: the 8 um positions are not on a regular grid")
    arr_max = (int(full.array_row.max()), int(full.array_col.max()))
    for b in sorted(set(bins) | {8}):
        f = b // inj.BASE_UM
        outd = os.path.join(views, f"square_{b:03d}um")
        if f == 1:
            m_out, meta = None, pos8[["barcode", "array_row", "array_col"]].copy()
            meta["n_children"] = 1
        else:
            m_out, meta = inj.aggregate(mat, pos8, f)
        A = np.column_stack([meta.array_row.to_numpy(float), meta.array_col.to_numpy(float), np.ones(len(meta))])
        xy = A @ coef
        meta["pxl_row_in_fullres"], meta["pxl_col_in_fullres"] = xy[:, 0], xy[:, 1]
        meta["in_tissue"] = 1
        meta["is_truncated"] = inj.truncated_mask(meta, arr_max[0], arr_max[1], f)
        os.makedirs(os.path.join(outd, "spatial"), exist_ok=True)
        if m_out is None:
            os.symlink(os.path.abspath(os.path.join(d8, "filtered_feature_bc_matrix")),
                       os.path.join(outd, "filtered_feature_bc_matrix"))
        else:
            inj.write_mex(os.path.join(outd, "filtered_feature_bc_matrix"), m_out, meta.barcode, feat)
        meta.to_parquet(os.path.join(outd, "spatial", "tissue_positions.parquet"), index=False)
        sf = dict(sf8); sf["bin_size_um"] = float(b)
        sf["spot_diameter_fullres"] = sf8["spot_diameter_fullres"] * f
        json.dump(sf, open(os.path.join(outd, "spatial", "scalefactors_json.json"), "w"), indent=4)
    return pos8, full

def metrics(mat, feat):
    mito = feat.symbol.astype(str).str.upper().str.startswith("MT-").to_numpy()
    tot = np.asarray(mat.sum(0)).ravel()
    det = np.diff(mat.tocsc().indptr)
    mt = np.asarray(mat[np.where(mito)[0]].sum(0)).ravel()
    frac = np.divide(mt, tot, out=np.zeros_like(tot, float), where=tot > 0)
    return tot, det, frac, float(tot.mean()), float(frac[tot > 0].mean()), int(mito.sum())

def boundary_band(full, bc8):
    from scipy.ndimage import distance_transform_edt
    tis = full.in_tissue.to_numpy() == 1
    g = np.zeros((int(full.array_row.max()) + 1, int(full.array_col.max()) + 1), bool)
    g[full.array_row.to_numpy()[tis], full.array_col.to_numpy()[tis]] = True
    d = distance_transform_edt(g)
    near = d[full.array_row.to_numpy(), full.array_col.to_numpy()]
    band = set(full.barcode.to_numpy()[tis & (near > 0) & (near <= EDGE_RADIUS_BINS)])
    return np.array([b in band for b in bc8])

def draw(out, pos8, tot, frac, removed_by, tests):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x, y = pos8.pxl_col_in_fullres.to_numpy(), pos8.pxl_row_in_fullres.to_numpy()
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.4))
    s = max(0.05, min(2.0, 2.5e5 / len(x)))
    for a, v, t, cm in ((ax[0], np.log10(tot + 1), "log10 UMIs per 8 um bin", "viridis"),
                        (ax[1], 100 * frac, "mitochondrial %", "magma")):
        sc = a.scatter(x, y, c=v, s=s, cmap=cm, linewidths=0, rasterized=True)
        fig.colorbar(sc, ax=a, shrink=0.7); a.set_title(t)
    any_rm = np.zeros(len(x), bool)
    ax[2].scatter(x, y, c="#d9d9d9", s=s, linewidths=0, rasterized=True)
    for (name, mask), col in zip(removed_by.items(), ("#d1495b", "#2e86ab", "#edae49")):
        ax[2].scatter(x[mask], y[mask], c=col, s=s, linewidths=0, rasterized=True,
                      label=f"{name}: {mask.sum():,} bins")
        any_rm |= mask
    ax[2].set_title(f"removed: {any_rm.mean():.1%} of in-tissue bins")
    ax[2].legend(loc="lower left", fontsize=8, markerscale=12, frameon=False)
    for a in ax:
        a.set_aspect("equal"); a.invert_yaxis(); a.set_xticks([]); a.set_yticks([])
    fig.suptitle("Image orientation: x and y are full-resolution image pixels", fontsize=9, y=0.02)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "qc_map.png"), dpi=150)

def main():
    ap = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Space Ranger outs/binned_outputs (holds square_008um)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--damage-um", type=float, default=None,
                    help="size of the damage seen on the image, in um; leave out if none is seen or its size is unknown")
    ap.add_argument("--edge", action="store_true", help="the damage sits at the tissue edge")
    ap.add_argument("--settings", default=None, help="guideline_settings.json (default: the one beside Fig. 6)")
    ap.add_argument("--keep-views", action="store_true", help="keep the coarse-bin views written for the tests")
    a = ap.parse_args()
    t0 = time.time()
    check_core()
    sp = settings_path(a); S = json.load(open(sp))
    d8 = os.path.join(a.input, "square_008um")
    if not os.path.exists(os.path.join(d8, "filtered_feature_bc_matrix", "matrix.mtx.gz")):
        sys.exit(f"{d8}/filtered_feature_bc_matrix/matrix.mtx.gz not found: --input is Space Ranger's binned_outputs")
    os.makedirs(a.out, exist_ok=True)
    views = os.path.join(a.out, "_views")
    if os.path.exists(views):
        shutil.rmtree(views)

    print("reading the 8 um bins", flush=True)
    mat, bc8, feat = inj.read_mex(os.path.join(d8, "filtered_feature_bc_matrix"))
    tot, det, frac, mean_umi, mean_mito, n_mito = metrics(mat, feat)
    tier, tests, notes = plan(S, mean_umi, a.damage_um, a.edge)
    if any(t["test"] == "spatialartifacts" for t in tests) and not os.path.exists(
            os.path.join(d8, "raw_feature_bc_matrix", "matrix.mtx.gz")):
        sys.exit("SpatialArtifacts needs square_008um/raw_feature_bc_matrix (the whole capture area)")
    for t in tests:
        print(f"  {t['test']} on {t['bin']} um bins" + (f", {t['k']} neighbours" if t["k"] else "") + f": {t['why']}")
    for n in notes:
        print("  note:", n)

    print("writing the bin views", flush=True)
    bins = {t["bin"] for t in tests if t["test"] != "spatialartifacts"}
    bins |= {16 for t in tests if t["test"] == "spatialartifacts" and t["bin"] == 16}
    pos8, full = prepare(a.input, views, bins, mat, bc8, feat)
    del mat

    removed_by = {}
    for t in tests:
        name = f"{t['test']}_{t['bin']}um" + (f"_k{t['k']}" if t["k"] else "")
        f = os.path.join(a.out, f"removed_{name}.txt")
        cmd = ["Rscript", os.path.join(HERE, "qc_run.R"), "--views", views, "--test", t["test"],
               "--bin", str(t["bin"]), "--out", f, "--source8", os.path.abspath(d8)]
        if t["k"]:
            cmd += ["--k", str(t["k"])]
        print("running", name, flush=True)
        subprocess.run(cmd, check=True)
        rm = set(l.strip() for l in open(f) if l.strip())
        removed_by[name] = np.array([b in rm for b in bc8])

    any_rm = np.zeros(len(bc8), bool)
    for m in removed_by.values():
        any_rm |= m
    band = boundary_band(full, bc8)
    with open(os.path.join(a.out, "removed_8um.txt"), "w") as fh:
        fh.write("".join(b + "\n" for b in np.asarray(bc8)[any_rm]))
    tab = pd.DataFrame({"barcode": bc8, "array_row": pos8.array_row, "array_col": pos8.array_col,
                        "sum": tot, "detected": det, "mito_frac": frac,
                        "truncated": inj.truncated_mask(pos8, int(full.array_row.max()), int(full.array_col.max()), 1),
                        "boundary_band": band, "removed": any_rm})
    for n, m in removed_by.items():
        tab[f"removed_{n}"] = m
    tab.to_csv(os.path.join(a.out, "bins_8um.csv.gz"), index=False)
    draw(a.out, pos8, tot, frac, removed_by, tests)

    rep = {"input": os.path.abspath(a.input), "settings": sp, "settings_generated": S.get("generated"),
           "damage_um": a.damage_um, "edge": a.edge, "in_tissue_bins": int(len(bc8)),
           "mean_umi_per_bin": round(mean_umi, 2), "mean_mito_fraction": round(mean_mito, 5),
           "mito_genes": n_mito, "depth_tier": tier, "tests": tests, "notes": notes,
           "removed_bins": int(any_rm.sum()), "removed_share": round(float(any_rm.mean()), 5),
           "removed_by_test": {n: int(m.sum()) for n, m in removed_by.items()},
           "boundary_band_share_of_tissue": round(float(band.mean()), 5),
           "boundary_band_share_of_removed": round(float(band[any_rm].mean()), 5) if any_rm.any() else None,
           "minutes": round((time.time() - t0) / 60, 1)}
    json.dump(rep, open(os.path.join(a.out, "qc_report.json"), "w"), indent=1)
    md = [f"# QC report", "", f"- input: `{rep['input']}`",
          f"- in-tissue 8 um bins: {len(bc8):,}",
          f"- baselines: {mean_umi:.1f} UMIs per bin, mean mitochondrial fraction {100 * mean_mito:.2f}% "
          f"({n_mito} mitochondrial genes)",
          f"- damage: " + (f"{a.damage_um:g} um" + (", at the tissue edge" if a.edge else "") if a.damage_um else "none given"),
          "", "## Tests (Fig. 6)", ""]
    md += [f"- **{t['test']}**, {t['bin']} um bins" + (" + A&D" if t["bin"] != 8 else "")
           + (f", {t['k']} neighbours" if t["k"] else "") + f": {t['why']}; removed {removed_by[n].sum():,} bins"
           for t, n in zip(tests, removed_by)]
    md += ["", "## What was removed", "",
           f"- {any_rm.sum():,} bins, {any_rm.mean():.1%} of the tissue",
           f"- the boundary band (within {EDGE_RADIUS_BINS * 8} um of the tissue edge) holds "
           f"{band.mean():.1%} of the tissue and "
           + (f"{band[any_rm].mean():.1%} of the removed bins" if any_rm.any() else "none of the removed bins"),
           "- check on qc_map.png that removed bins sit on the damage, not on one cell type or along the edge",
           ] + [f"- note: {n}" for n in notes]
    open(os.path.join(a.out, "qc_report.md"), "w").write("\n".join(md) + "\n")
    if not a.keep_views:
        shutil.rmtree(views)
    print(f"done in {rep['minutes']} min: {any_rm.sum():,} of {len(bc8):,} bins removed ({any_rm.mean():.1%})")

if __name__ == "__main__":
    main()
