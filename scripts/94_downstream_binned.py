import argparse, os, re
import numpy as np, pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from scipy.sparse import csr_matrix, coo_matrix
from sklearn.metrics import adjusted_rand_score

sc.settings.verbosity = 0
BC = re.compile(r"s_(\d+)um_(\d+)_(\d+)")

def read_bins(cell):
    d = os.path.join(cell, "clean", "square_008um", "filtered_feature_bc_matrix")
    m = mmread(os.path.join(d, "matrix.mtx.gz")).tocsc().T.tocsr()
    bc = pd.read_csv(os.path.join(d, "barcodes.tsv.gz"), header=None)[0].to_numpy()
    ft = pd.read_csv(os.path.join(d, "features.tsv.gz"), sep="\t", header=None)
    a = ad.AnnData(csr_matrix(m))
    a.obs_names = pd.Index(bc, name=None)
    a.var_names = pd.Index(ft[1].astype(str).to_numpy(), name=None)
    a.var_names_make_unique()
    return a

def coarse_labels(barcodes, bin_um):
    step = bin_um // 8
    out = np.empty(len(barcodes), dtype=object)
    for i, b in enumerate(barcodes):
        m = BC.match(b)
        r, c = int(m.group(2)) // step, int(m.group(3)) // step
        out[i] = f"r{r}_c{c}"
    return out

def pool(a, bin_um):
    lab = coarse_labels(a.obs_names.to_numpy(), bin_um)
    uniq, inv = np.unique(lab, return_inverse=True)
    G = coo_matrix((np.ones(len(inv)), (inv, np.arange(len(inv)))),
                   shape=(len(uniq), a.n_obs)).tocsr()
    X = G @ a.X
    out = ad.AnnData(csr_matrix(X))
    out.obs_names = pd.Index(uniq, name=None)
    out.var_names = pd.Index(a.var_names.to_numpy(), name=None)
    keep = np.asarray(out.X.sum(axis=1)).ravel() > 0
    return out[keep].copy()

def cluster(a, seed, npc=30, k=15, res=1.0):
    a = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    a = a[:, a.var.highly_variable].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=npc, svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(a, n_neighbors=k, n_pcs=npc, random_state=seed)
    sc.tl.leiden(a, resolution=res, key_added="cl", random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    return pd.Series(a.obs["cl"].astype(str).to_numpy(), index=a.obs_names)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="/work/benchmark_grid")
    ap.add_argument("--section", required=True)
    ap.add_argument("--deep", type=int, default=100)
    ap.add_argument("--shallow", type=int, default=12)
    ap.add_argument("--bin", type=int, default=32, help="clustering bin width in um")
    ap.add_argument("--arms", default="01_no_qc,02_mad_8um,10_mad_64um,03_spotsweeper_8um")
    ap.add_argument("--out", default="/work/benchmark_grid/figures/downstream_binned.csv")
    a = ap.parse_args()

    rows = []
    cells = {d: os.path.join(a.grid, f"{a.section}_d{d}") for d in (a.deep, a.shallow)}
    raw = {d: read_bins(c) for d, c in cells.items()}
    for d, x in raw.items():
        print(f"{a.section} d{d}: {x.n_obs} bins of 8 µm, {x.n_vars} genes", flush=True)

    for arm in a.arms.split(","):
        lab = {}
        for dep, cell in cells.items():
            f = os.path.join(cell, "results", arm, "clean.txt")
            drop = set(open(f).read().split()) if os.path.exists(f) else set()
            keep = [b for b in raw[dep].obs_names if b not in drop]
            sub = pool(raw[dep][keep].copy(), a.bin)
            if sub.n_obs < 500:
                print(f"  {arm} d{dep}: only {sub.n_obs} bins, skipped", flush=True)
                lab = {}; break
            lab[dep] = cluster(sub, seed=0)
            if dep == a.deep:
                reseed = cluster(sub, seed=1)
                sh = lab[dep].index.intersection(reseed.index)
                rows.append(dict(section=a.section, arm=arm, comparison="reseed",
                                 bin_um=a.bin, depth_a=dep, depth_b=dep, n_shared=len(sh),
                                 ari=adjusted_rand_score(lab[dep][sh], reseed[sh]),
                                 n_clusters_a=lab[dep].nunique(),
                                 n_clusters_b=reseed.nunique(),
                                 kept_a=sub.n_obs, kept_b=sub.n_obs))
                print(f"  {arm} d{dep}: {len(keep)} of {raw[dep].n_obs} 8 µm bins kept, "
                      f"{sub.n_obs} bins of {a.bin} µm, {lab[dep].nunique()} clusters, "
                      f"reseed ARI {rows[-1]['ari']:.3f}", flush=True)
            else:
                print(f"  {arm} d{dep}: {len(keep)} of {raw[dep].n_obs} 8 µm bins kept, "
                      f"{sub.n_obs} bins of {a.bin} µm, {lab[dep].nunique()} clusters",
                      flush=True)
        if len(lab) != 2:
            continue
        sh = lab[a.deep].index.intersection(lab[a.shallow].index)
        rows.append(dict(section=a.section, arm=arm, comparison="across_depth",
                         bin_um=a.bin, depth_a=a.deep, depth_b=a.shallow, n_shared=len(sh),
                         ari=adjusted_rand_score(lab[a.deep][sh], lab[a.shallow][sh]),
                         n_clusters_a=lab[a.deep].nunique(),
                         n_clusters_b=lab[a.shallow].nunique(),
                         kept_a=len(lab[a.deep]), kept_b=len(lab[a.shallow])))
        print(f"  -> {arm} ARI across depth {rows[-1]['ari']:.3f} "
              f"on {len(sh)} shared bins", flush=True)

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    df.to_csv(a.out, mode="a", header=not os.path.exists(a.out), index=False)
    print("appended", len(df), "rows to", a.out)

if __name__ == "__main__":
    main()
