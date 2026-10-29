import sys, os, gzip
out = []
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp
for run in sys.argv[1:]:
    d = f"/work/benchmark_grid/{run}/mito_64um"
    m = os.path.join(d, "square_008um", "filtered_feature_bc_matrix")
    X = scipy.io.mmread(gzip.open(os.path.join(m, "matrix.mtx.gz"))).tocsc()
    feat = pd.read_csv(os.path.join(m, "features.tsv.gz"), sep="\t", header=None)
    bc = pd.read_csv(os.path.join(m, "barcodes.tsv.gz"), header=None)[0]
    mt = feat[1].astype(str).str.upper().str.startswith("MT-").to_numpy()
    tot = np.asarray(X.sum(0)).ravel(); mtc = np.asarray(X[mt].sum(0)).ravel()
    frac = np.divide(mtc, tot, out=np.zeros_like(tot, float), where=tot > 0) * 100
    g = pd.read_parquet(os.path.join(d, "ground_truth_008um.parquet"))
    col = "barcode" if "barcode" in g.columns else g.columns[0]
    gt = set(g.barcode[g.is_artifact])
    inside = bc.isin(gt).to_numpy()
    med = np.median(frac); mad = 1.4826 * np.median(np.abs(frac - med))
    out.append(f"{run}: bins {len(bc)}, LQR bins {inside.sum()}; mito % outside median {np.median(frac[~inside]):.1f} "
          f"(p90 {np.percentile(frac[~inside], 90):.1f}), inside median {np.median(frac[inside]):.1f}; "
          f"3-MAD upper threshold {med + 3 * mad:.1f}% (median {med:.1f}, scaled MAD {mad:.1f}); "
          f"inside above threshold {np.mean(frac[inside] > med + 3 * mad):.1%}")

open("/work/benchmark_grid/figures/mito_threshold.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
