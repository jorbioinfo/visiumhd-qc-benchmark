#!/usr/bin/env python3
import gzip, os, sys
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import binomtest

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(WK, "scripts"))
import panel as P

GRID = os.path.join(WK, "benchmark_grid")
OUT = os.path.join(GRID, "figures")
K, CUT, SS_ARM = 48, 3.0, "24_spotsweeper_8um_k48"

def read_mtx(base, genes=False):
    with gzip.open(f"{base}/matrix.mtx.gz", "rt") as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        ng, nb, nnz = map(int, line.split())
        m = pd.read_csv(fh, sep=" ", header=None, names=["g", "b", "v"],
                        dtype={"g": np.int32, "b": np.int32, "v": np.int32}, engine="c")
    feat = pd.read_csv(f"{base}/features.tsv.gz", sep="\t", header=None)
    mito_g = np.where(feat[1].astype(str).str.match(r"(?i)^MT-"))[0] + 1
    b = m.b.values - 1
    s = np.bincount(b, weights=m.v.values, minlength=nb)
    det = np.bincount(b, minlength=nb).astype(float)
    mm = np.isin(m.g.values, mito_g)
    mt = np.bincount(b[mm], weights=m.v.values[mm], minlength=nb)
    bc = pd.read_csv(f"{base}/barcodes.tsv.gz", header=None)[0].values
    return bc, s, det, mt

def pct_of(s, mt):
    return np.where(s > 0, 100 * mt / np.where(s > 0, s, 1), 0.0)

def lower(x):
    with np.errstate(divide="ignore", invalid="ignore"):
        y = np.log2(x)
        med = np.median(y)
        mad = 1.4826 * np.median(np.abs(y - med))
        return y < med - 3 * mad

def higher(x):
    med = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - med))
    return x > med + 3 * mad, med, mad

def modz(v):
    med = np.median(v, axis=1, keepdims=True)
    mad = 1.4826 * np.median(np.abs(v - med), axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = 0.6745 * (v[:, :1] - med) / mad
    z = z[:, 0].copy()
    z[~np.isfinite(z)] = 0.0
    return z, mad[:, 0], v[:, 0] - med[:, 0]

def removed(cell, arm, cond):
    with open(f"{cell}/results/{arm}/{cond}.txt") as fh:
        return {l.strip() for l in fh if l.strip()}

def check_mad():
    rows = []
    for sec, dep, cell in P.cells(GRID, (100, 12)):
        base = f"{cell}/clean/square_008um"
        bc, s, det, mt = read_mtx(f"{base}/filtered_feature_bc_matrix")
        pos = pd.read_parquet(f"{base}/spatial/tissue_positions.parquet").set_index("barcode")
        keep = pos.loc[bc, "in_tissue"].values == 1
        s, det, mt = s[keep], det[keep], mt[keep]
        pct = pct_of(s, mt)
        lo_s, lo_d = lower(s), lower(det)
        hi_m, mmed, mmad = higher(pct)
        union = lo_s | lo_d | hi_m
        sc = pd.read_csv(f"{cell}/scored/scores.csv")
        arm = sc[(sc.arm == "02_mad_8um") & (sc.condition == "clean")].discard_frac.iloc[0]
        rows.append(dict(section=P.short(sec), chem=P.chemistry(sec), depth=dep, n=int(keep.sum()),
                         arm_fpr=100 * arm, repl_fpr=100 * union.mean(), empty=100 * (s == 0).mean(),
                         libdet=100 * (lo_s | lo_d).mean(), mito_only=100 * (hi_m & ~lo_s & ~lo_d).mean(),
                         mito_med=mmed, mito_mad=mmad, mito_zero=100 * (mt == 0).mean()))
        print("mad", sec, dep, flush=True)
    return pd.DataFrame(rows)

def check_spot():
    rows = []
    for sec, dep, cell in P.cells(GRID, (100, 12)):
        for cond in ("mito_32um", "dropout01_32um"):
            base = f"{cell}/{cond}/square_008um"
            bc, s, det, mt = read_mtx(f"{base}/filtered_feature_bc_matrix")
            pos = pd.read_parquet(f"{base}/spatial/tissue_positions.parquet").set_index("barcode").loc[bc]
            keep = ~pos.is_truncated.values.astype(bool)
            bc, s, det, mt = bc[keep], s[keep], det[keep], mt[keep]
            xy = pos[["array_row", "array_col"]].values[keep].astype(float)
            pct = pct_of(s, mt)
            gt = pd.read_parquet(f"{cell}/{cond}/ground_truth_008um.parquet")
            li = np.where(np.isin(bc, gt.barcode[gt.is_artifact].values))[0]
            actual = np.isin(bc[li], list(removed(cell, SS_ARM, cond)))
            _, nn = cKDTree(xy).query(xy[li], k=K + 1)
            zs, ms, ds = modz(np.log1p(s)[nn])
            zd, _, _ = modz(np.log1p(det)[nn])
            zm, mm_, _ = modz(pct[nn])
            pred = (zs < -CUT) | (zd < -CUT) | (zm > CUT)
            miss = ~actual
            r = dict(section=P.short(sec), chem=P.chemistry(sec), depth=dep, cond=cond, n_lqr=len(li),
                     recall=100 * actual.mean(), agree=100 * (pred == actual).mean())
            if cond.startswith("mito"):
                r.update(mito_mad0=100 * (mm_ == 0).mean(),
                         missed_mad0=100 * ((mm_ == 0) & miss).sum() / max(miss.sum(), 1))
            else:
                r.update(lib_gap=float(np.median(ds)), lib_mad=float(np.median(ms)))
            rows.append(r)
        print("spot", sec, dep, flush=True)
    return pd.DataFrame(rows)

def check_ap():
    rows = []
    for sec, dep, cell in P.cells(GRID, (12,)):
        r = dict(section=P.short(sec), chem=P.chemistry(sec), depth=dep)
        for size in (8, 16, 32, 64):
            base = f"{cell}/clean/square_{size:03d}um"
            bc, s, _, mt = read_mtx(f"{base}/filtered_feature_bc_matrix")
            pos = pd.read_parquet(f"{base}/spatial/tissue_positions.parquet").set_index("barcode").loc[bc]
            keep = (pos.in_tissue.values == 1) & ~pos.is_truncated.values.astype(bool)
            pct = pct_of(s[keep], mt[keep])
            r[f"mito_zero{size}"] = 100 * float((mt[keep] == 0).mean())
            r[f"mito_med{size}"] = float(np.median(pct))
            r[f"mito_mad{size}"] = float(np.median(np.abs(pct - np.median(pct))))
            r[f"empty{size}"] = 100 * float((s[keep] == 0).mean())
        rows.append(r)
        print("ap", sec, flush=True)
    return pd.DataFrame(rows)

def sign(x, y):
    d = (x - y).dropna()
    up, dn = int((d > 0).sum()), int((d < 0).sum())
    return up, dn, (binomtest(up, up + dn).pvalue if up + dn else float("nan"))

def report(M, S, A):
    print("\n## mad: the MAD threshold on 8 um bins, uninjected condition")
    print(f"recomputed against the arm: largest difference {np.abs(M.repl_fpr - M.arm_fpr).max():.2f} points")
    for ch in ("probe", "3prime"):
        for dep in (100, 12):
            E = M[(M.chem == ch) & (M.depth == dep)]
            z = ((E.mito_med == 0) & (E.mito_mad == 0)).sum()
            print(f"{ch:6} d{dep:<3} n={len(E):2} fpr {E.arm_fpr.median():.1f}%: mito only {E.mito_only.median():.1f}, "
                  f"library size or genes {E.libdet.median():.1f} (empty bins {E['empty'].median():.1f}); "
                  f"bins without a mito molecule {E.mito_zero.median():.1f}%; median and MAD of mito % zero in {z} of {len(E)}")
        W = M[M.chem == ch].pivot(index="section", columns="depth")
        rise = W.arm_fpr[12] - W.arm_fpr[100]
        dm, dl = W.mito_only[12] - W.mito_only[100], W.libdet[12] - W.libdet[100]
        up, dn, p = sign(dm, dl)
        share = (dm / (dm + dl)).where(rise > 0)
        print(f"{ch:6} rise 100 -> 12: {rise.median():.1f} points, mito only {dm.median():.1f}, library size or genes "
              f"{dl.median():.1f}; mito part rose more in {up} of {up + dn} (p = {p:.2g}); share of the rise from mito "
              f"only, median {100 * share.median():.0f}%")
    print("\n## spot: SpotSweeper, 48 neighbours, 8 um bins, 32 um-diameter LQRs")
    print(f"recomputed against the arm: verdict agreement on LQR bins, lowest {S.agree.min():.1f}%, median {S.agree.median():.1f}%")
    for ch in ("probe", "3prime"):
        for dep in (100, 12):
            E = S[(S.chem == ch) & (S.depth == dep)]
            mi, sl = E[E.cond == "mito_32um"], E[E.cond == "dropout01_32um"]
            print(f"{ch:6} d{dep:<3} n={len(mi):2} mito: recall {mi.recall.median():.1f}%, LQR bins with neighbourhood MAD 0 "
                  f"{mi.mito_mad0.median():.1f}%, of missed bins {mi.missed_mad0.median():.1f}%; signal loss: recall "
                  f"{sl.recall.median():.1f}%, gap to neighbourhood median {sl.lib_gap.median():.2f} log units, "
                  f"neighbourhood MAD {sl.lib_mad.median():.2f}")
        mi = S[(S.chem == ch) & (S.cond == "mito_32um")].pivot(index="section", columns="depth", values="mito_mad0")
        up, dn, p = sign(mi[12], mi[100])
        print(f"{ch:6} mito LQR bins with MAD 0, 12 against 100: higher in {up} of {up + dn} (p = {p:.2g})")
        sl = S[(S.chem == ch) & (S.cond == "dropout01_32um")].pivot(index="section", columns="depth")
        up, dn, p = sign(sl.lib_gap[12], sl.lib_gap[100])
        print(f"{ch:6} signal-loss gap, 12 against 100: closer to zero in {up} of {up + dn} (p = {p:.2g})")
        up, dn, p = sign(sl.lib_mad[12], sl.lib_mad[100])
        print(f"{ch:6} signal-loss neighbourhood MAD, 12 against 100: wider in {up} of {up + dn} (p = {p:.2g})")
    print("\n## ap: 12 UMIs per bin, uninjected condition")
    for ch in ("probe", "3prime"):
        E = A[A.chem == ch]
        for size in (8, 16, 32, 64):
            z = ((E[f"mito_med{size}"] == 0) & (E[f"mito_mad{size}"] == 0)).sum()
            print(f"{ch:6} n={len(E):2} {size:2} um: bins without a mito molecule {E[f'mito_zero{size}'].median():.1f}%, "
                  f"median and MAD of mito % zero in {z}, empty bins {E[f'empty{size}'].median():.2f}%")

if __name__ == "__main__":
    which = sys.argv[1:] or ["mad", "spot", "ap"]
    T = {}
    for k, f in (("mad", check_mad), ("spot", check_spot), ("ap", check_ap)):
        path = os.path.join(OUT, f"native_checks_{k}.csv")
        if k in which:
            f().to_csv(path, index=False)
        T[k] = pd.read_csv(path)
    report(T["mad"], T["spot"], T["ap"])
