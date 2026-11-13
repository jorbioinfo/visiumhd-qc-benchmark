#!/usr/bin/env python3
import argparse, importlib.util, json, os, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

spec = importlib.util.spec_from_file_location("inj", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                  "10_inject_artifacts.py"))
inj = importlib.util.module_from_spec(spec); spec.loader.exec_module(inj)
_gzip_open = inj.gzip.open
inj.gzip.open = lambda *a, **k: _gzip_open(*a, compresslevel=1, **k)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--depths", type=float, nargs="+", required=True)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    d8 = f"{a.src}/square_008um"
    mat0, bc8, feat = inj.read_mex(f"{d8}/filtered_feature_bc_matrix")
    pos8 = pq.read_table(f"{d8}/spatial/tissue_positions.parquet").to_pandas()
    sf8 = json.load(open(f"{d8}/spatial/scalefactors_json.json"))
    native = mat0.sum() / mat0.shape[1]
    pos8 = pos8.set_index("barcode").reindex(bc8).reset_index()
    assert pos8.barcode.tolist() == list(bc8)
    full_pos8 = pq.read_table(f"{d8}/spatial/tissue_positions.parquet").to_pandas()
    coef, resid = inj.grid_basis(full_pos8)
    arr_max = (int(full_pos8.array_row.max()), int(full_pos8.array_col.max()))
    inj.log(f"{a.name}: native {native:.1f} UMIs/bin, {mat0.shape[1]:,} bins, grid residual {resid:.4f} px")

    for t in a.depths:
        tag = f"{a.name}_d{t:g}"
        out = os.path.join(a.out_root, tag)
        if t > native:
            inj.log(f"skip {tag}: target exceeds native depth {native:.1f}")
            continue
        if os.path.exists(os.path.join(out, "manifest.csv")):
            inj.log(f"skip {tag}: already built")
            continue
        rng = np.random.default_rng(a.seed)
        p = t / native
        mat = mat0.tocsc(copy=True)
        if p < 1.0:
            mat.data = rng.binomial(mat.data.astype(np.int64), p).astype(mat.dtype)
            mat.eliminate_zeros()
        inj.log(f"{tag}: fraction {p:.4f}, {mat.sum() / mat.shape[1]:.1f} UMIs/bin")
        root = os.path.join(out, "clean")
        os.makedirs(root, exist_ok=True)
        pd.DataFrame({"barcode": bc8, "is_artifact": np.zeros(len(bc8), bool)}).to_parquet(
            os.path.join(root, "ground_truth_008um.parquet"), index=False)
        manifest = []
        for b in inj.BIN_SIZES:
            f = b // inj.BASE_UM
            outd = os.path.join(root, f"square_{b:03d}um")
            if f == 1:
                m_out, meta = mat, pos8[["barcode", "array_row", "array_col"]].copy()
                meta["n_children"] = 1
            else:
                m_out, meta = inj.aggregate(mat, pos8, f)
            A = np.column_stack([meta.array_row.to_numpy(float), meta.array_col.to_numpy(float), np.ones(len(meta))])
            xy = A @ coef
            meta["pxl_row_in_fullres"] = xy[:, 0]
            meta["pxl_col_in_fullres"] = xy[:, 1]
            meta["in_tissue"] = 1
            meta["is_truncated"] = inj.truncated_mask(meta, arr_max[0], arr_max[1], f)
            inj.write_mex(os.path.join(outd, "filtered_feature_bc_matrix"), m_out, meta.barcode, feat)
            os.makedirs(os.path.join(outd, "spatial"), exist_ok=True)
            meta.to_parquet(os.path.join(outd, "spatial", "tissue_positions.parquet"), index=False)
            sf = dict(sf8); sf["bin_size_um"] = float(b)
            sf["spot_diameter_fullres"] = sf8["spot_diameter_fullres"] * f
            json.dump(sf, open(os.path.join(outd, "spatial", "scalefactors_json.json"), "w"), indent=4)
            manifest.append({"condition": "clean", "kind": "clean", "diameter_um": None, "bin_um": b,
                             "n_bins": int(m_out.shape[1]), "n_truncated": int(meta.is_truncated.sum()),
                             "path": outd})
        pd.DataFrame(manifest).to_csv(os.path.join(out, "manifest.csv"), index=False)
        json.dump(dict(seed=a.seed, target_umi_per_bin=t, subsample=p, native_umi_per_bin=float(native),
                       src=a.src, built_by="145_molecule_curve_build.py"),
                  open(os.path.join(out, "injection_params.json"), "w"), indent=2)

if __name__ == "__main__":
    main()
