import argparse, gzip, os, sys
import numpy as np, pandas as pd, scipy.io
import pyarrow.parquet as pq
from scipy import stats

N_PANEL = 20
TILE = 32
MIN_PER_GROUP = 10

def load(cond_dir):
    d = os.path.join(cond_dir, "square_008um", "filtered_feature_bc_matrix")
    m = scipy.io.mmread(os.path.join(d, "matrix.mtx.gz")).tocsr()
    bc = np.array([l.strip() for l in gzip.open(os.path.join(d, "barcodes.tsv.gz"), "rt")])
    ft = pd.read_csv(os.path.join(d, "features.tsv.gz"), sep="\t", header=None,
                     names=["id", "symbol", "type"])
    pos = pq.read_table(os.path.join(cond_dir, "square_008um", "spatial",
                                     "tissue_positions.parquet")).to_pandas()
    pos = pos.set_index("barcode").reindex(bc)
    return m, bc, ft, pos

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", required=True)
    ap.add_argument("--name", required=True, help="label written into the output")
    ap.add_argument("--n-panel", type=int, default=N_PANEL)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    m, bc, ft, pos = load(os.path.join(a.benchmark, "clean"))
    mcsc = m.tocsc()
    total = np.asarray(mcsc.sum(0)).ravel()
    is_mt = ft.symbol.astype(str).str.upper().str.startswith("MT-").to_numpy()
    order = np.argsort(-np.asarray(mcsc.sum(1)).ravel())
    rows_idx = np.array([i for i in order if not is_mt[i]][:a.n_panel])
    mark = np.asarray(mcsc[rows_idx, :].sum(0)).ravel()
    panel_share = mark.sum() / max(total.sum(), 1)
    print(f"{a.name}: panel = top {a.n_panel} non-mito genes, {panel_share:.1%} of counts")
    print("  " + ", ".join(ft.symbol.iloc[rows_idx[:10]].astype(str)))

    keep = total > 0
    frac = np.zeros(len(bc)); frac[keep] = mark[keep] / total[keep]
    r, c = pos.array_row.to_numpy(), pos.array_col.to_numpy()
    tile = (r // TILE).astype(np.int64) * 100000 + (c // TILE)
    lg = np.log10(total + 1)
    edges = np.percentile(lg[keep], np.linspace(0, 100, 21))
    dec = np.clip(np.digitize(lg, edges[1:-1]), 0, 19)

    out = []
    res = os.path.join(a.benchmark, "results")
    for arm in sorted(os.listdir(res)):
        f = os.path.join(res, arm, "clean.txt")
        if arm == "logs" or not os.path.isfile(f):
            continue
        flagged = set(l.strip() for l in open(f))
        if not flagged:
            continue
        fl = np.isin(bc, list(flagged)) & keep
        un = (~fl) & keep
        if fl.sum() < 50 or un.sum() < 50:
            continue

        naive = np.median(frac[fl]) - np.median(frac[un])

        diffs, wts = [], []
        for dv in np.unique(dec[keep]):
            s = (dec == dv)
            if (fl & s).sum() >= 20 and (un & s).sum() >= 20:
                diffs.append(np.median(frac[fl & s]) - np.median(frac[un & s]))
                wts.append((fl & s).sum())
        by_depth = np.average(diffs, weights=wts) if diffs else np.nan

        by_tile, p_tile, ratio, n_cells = np.nan, np.nan, np.nan, 0
        scale_um, n_strata = np.nan, np.nan
        for side, nstr in ((TILE, 20), (TILE * 2, 10), (TILE * 4, 5)):
            tl = (r // side).astype(np.int64) * 100000 + (c // side)
            if nstr == 20:
                dc = dec
            else:
                e = np.percentile(lg[keep], np.linspace(0, 100, nstr + 1))
                dc = np.clip(np.digitize(lg, e[1:-1]), 0, nstr - 1)
            pf, pu = [], []
            for tv in np.unique(tl):
                st = (tl == tv)
                for dv in np.unique(dc[st & keep]):
                    s = st & (dc == dv)
                    if (fl & s).sum() >= MIN_PER_GROUP and (un & s).sum() >= MIN_PER_GROUP:
                        pf.append(np.median(frac[fl & s])); pu.append(np.median(frac[un & s]))
            if len(pf) < 5:
                continue
            pf, pu = np.array(pf), np.array(pu)
            by_tile = np.mean(pf - pu)
            try:
                p_tile = stats.wilcoxon(pf, pu).pvalue
            except ValueError:
                p_tile = np.nan
            ratio = np.mean(pf) / np.mean(pu) if np.mean(pu) > 0 else np.nan
            n_cells, scale_um, n_strata = len(pf), side * 8, nstr
            break

        out.append(dict(dataset=a.name, arm=arm, n_flagged=int(fl.sum()),
                        pct_flagged=100 * fl.sum() / keep.sum(),
                        marker_frac_flagged=np.median(frac[fl]),
                        marker_frac_unflagged=np.median(frac[un]),
                        diff_naive=naive, diff_by_depth=by_depth,
                        diff_matched=by_tile, ratio_matched=ratio,
                        n_cells=n_cells, match_scale_um=scale_um,
                        match_strata=n_strata, panel_share=panel_share,
                        p_tile=p_tile))
    df = pd.DataFrame(out)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    df.to_csv(a.out, index=False)
    print(df.to_string(index=False))
    print("wrote", a.out)

if __name__ == "__main__":
    main()
