import argparse, importlib.util, json, os, sys
import numpy as np, pandas as pd, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "inj", os.path.join(HERE, "10_inject_artifacts.py"))
inj = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inj)

def add_bin(cond_dir, bin_um, force=False):
    f = bin_um // inj.BASE_UM
    outd = os.path.join(cond_dir, f"square_{bin_um:03d}um")
    done = os.path.join(outd, "filtered_feature_bc_matrix", "matrix.mtx.gz")
    if os.path.exists(done) and not force:
        return f"{os.path.basename(cond_dir)}: {bin_um} um already present"
    d8 = os.path.join(cond_dir, "square_008um")
    mat, bc8, feat = inj.read_mex(os.path.join(d8, "filtered_feature_bc_matrix"))
    pos8 = pd.read_parquet(os.path.join(d8, "spatial", "tissue_positions.parquet"))
    pos8 = pos8.set_index("barcode").loc[bc8].reset_index()

    m_out, meta = inj.aggregate(mat, pos8, f)
    src = json.load(open(os.path.join(os.path.dirname(cond_dir), "injection_params.json")))["src"]
    full = pd.read_parquet(os.path.join(src, "square_008um", "spatial", "tissue_positions.parquet"))
    coef, resid = inj.grid_basis(full)
    if resid > 1.0:
        raise RuntimeError(f"grid basis residual {resid:.3f} px in {cond_dir}; "
                           "the 8 um positions are not on a regular grid")
    A = np.column_stack([meta.array_row.to_numpy(float),
                         meta.array_col.to_numpy(float), np.ones(len(meta))])
    xy = A @ coef
    meta["pxl_row_in_fullres"] = xy[:, 0]
    meta["pxl_col_in_fullres"] = xy[:, 1]
    meta["in_tissue"] = 1
    arr_max = (int(full.array_row.max()), int(full.array_col.max()))
    meta["is_truncated"] = inj.truncated_mask(meta, arr_max[0], arr_max[1], f)

    inj.write_mex(os.path.join(outd, "filtered_feature_bc_matrix"), m_out, meta.barcode, feat)
    os.makedirs(os.path.join(outd, "spatial"), exist_ok=True)
    meta.to_parquet(os.path.join(outd, "spatial", "tissue_positions.parquet"), index=False)
    sf8 = json.load(open(os.path.join(d8, "spatial", "scalefactors_json.json")))
    sf = dict(sf8)
    sf["bin_size_um"] = float(bin_um)
    sf["spot_diameter_fullres"] = sf8["spot_diameter_fullres"] * f
    json.dump(sf, open(os.path.join(outd, "spatial", "scalefactors_json.json"), "w"),
              indent=4)
    return (f"{os.path.basename(cond_dir)}: wrote {bin_um} um, {m_out.shape[1]:,} bins "
            f"({int(meta.is_truncated.sum())} truncated), "
            f"{int((meta.n_children < f * f).sum()):,} partially covered")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", required=True, help="a grid cell directory")
    ap.add_argument("--bin-um", type=int, default=128)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if a.bin_um % inj.BASE_UM or (a.bin_um // inj.BASE_UM) & ((a.bin_um // inj.BASE_UM) - 1):
        sys.exit(f"--bin-um must be 8 times a power of two, got {a.bin_um}")
    conds = [d for d in sorted(os.listdir(a.benchmark))
             if d not in ("results", "scored")
             and os.path.isdir(os.path.join(a.benchmark, d, "square_008um"))]
    if not conds:
        sys.exit(f"no conditions with square_008um under {a.benchmark}")
    for c in conds:
        print(add_bin(os.path.join(a.benchmark, c), a.bin_um, a.force), flush=True)

if __name__ == "__main__":
    main()
