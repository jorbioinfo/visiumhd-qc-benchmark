#!/usr/bin/env python3
import importlib.util, math, os, sys
import numpy as np, pandas as pd

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(WK, "scripts"))
import panel as P

GRID = os.path.join(WK, "benchmark_grid")
SAS = os.path.join(WK, "benchmark_saswitch")
spec = importlib.util.spec_from_file_location("m79", os.path.join(WK, "scripts", "79_build_model_tables.py"))
m79 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m79)

WIDE = {32: "44_spotsweeper_8um_k192", 64: "45_spotsweeper_8um_k768"}
AP = {32: "25_spotsweeper_16um_k48", 64: "26_spotsweeper_32um_k48"}
NATIVE = "24_spotsweeper_8um_k48"
SAME = {100: (8, "02_mad_8um", "19_spatialartifacts_8um"), 12: (16, "08_mad_16um", "20_spatialartifacts_16um")}
KINDS = ("dropout01edge", "mitoedge")

def sign_p(up, dn):
    n = up + dn
    if n == 0:
        return float("nan")
    return min(1.0, 2 * sum(math.comb(n, k) for k in range(min(up, dn) + 1)) / 2 ** n)

def test(x, y):
    i = x.dropna().index.intersection(y.dropna().index)
    x, y = x[i], y[i]
    up, dn = int((x > y).sum()), int((x < y).sum())
    return dict(n=len(i), a=100 * x.median(), b=100 * y.median(), higher=up, lower=dn, p=sign_p(up, dn))

def say(tag, r):
    print(f"  {tag:58s} {r['a']:5.1f} vs {r['b']:5.1f}  higher {r['higher']:2d} lower {r['lower']:2d} "
          f"of {r['n']:2d}  p {r['p']:.2g}")

def removed(path):
    if not os.path.exists(path):
        return None
    return {ln.strip() for ln in open(path) if ln.strip()}

def truth(cell, cond):
    d = os.path.join(cell, cond)
    gt = pd.read_parquet(os.path.join(d, "ground_truth_008um.parquet"), columns=["barcode", "is_artifact"])
    p = pd.read_parquet(os.path.join(d, "square_008um", "spatial", "tissue_positions.parquet"))
    uni = set(p.barcode[p.in_tissue == 1])
    return set(gt.barcode[gt.is_artifact]) & uni, uni

rows = []
S = []
for sec, dep, cell in P.cells(GRID, depths=(100, 12)):
    f = P.scores(cell)
    f["section"], f["depth"], f["chem"] = sec, dep, P.chemistry(sec)
    S.append(f)
S = pd.concat(S, ignore_index=True)
S["kind"] = S["kind"].astype(str)

def rec(arm, kind, dm, dep, chem):
    s = S[(S.arm == arm) & (S.kind == kind) & (S.diameter_um == dm) & (S.depth == dep) & (S.chem == chem)]
    return s.set_index("section").recall

print("A. wider neighbourhood against A&P at the edge, 100 UMIs per bin")
for chem in ("probe", "3prime"):
    print(f" {chem}")
    for dm in (32, 64):
        mean_w, mean_a = [], []
        for kind in KINDS:
            w, a, n8 = (rec(x, kind, dm, 100, chem) for x in (WIDE[dm], AP[dm], NATIVE))
            r = test(w, a); say(f"{kind} {dm} um: {WIDE[dm]} vs {AP[dm]}", r)
            rows.append(dict(part="A", chem=chem, depth=100, cond=f"{kind}_{dm}um", arm_a=WIDE[dm], arm_b=AP[dm], **r))
            r = test(w, n8); say(f"{kind} {dm} um: {WIDE[dm]} vs {NATIVE}", r)
            rows.append(dict(part="A", chem=chem, depth=100, cond=f"{kind}_{dm}um", arm_a=WIDE[dm], arm_b=NATIVE, **r))
            mean_w.append(w); mean_a.append(a)
        w = pd.concat(mean_w, axis=1).mean(axis=1, skipna=False)
        a = pd.concat(mean_a, axis=1).mean(axis=1, skipna=False)
        r = test(w, a); say(f"mean of the two types {dm} um: wider vs A&P", r)
        rows.append(dict(part="A", chem=chem, depth=100, cond=f"mean_{dm}um", arm_a=WIDE[dm], arm_b=AP[dm], **r))
for dm in (32, 64):
    for kind in KINDS:
        n = rec(WIDE[dm], kind, dm, 12, "probe").dropna().size
        print(f"  12 UMIs per bin, {kind} {dm} um, {WIDE[dm]}: {n} of 33 probe-based sections scored so far")

print("\nB, C, D from the removed barcodes")
U, M = [], []
for sec, dep, cell in P.cells(GRID, depths=(100, 12)):
    run, chem = os.path.basename(cell), P.chemistry(sec)
    b, mad, sa = SAME[dep]
    band = m79.edge_barcodes(GRID, run, 16.0)
    for cond in ["clean"] + [f"{k}_{dm}um" for k in KINDS for dm in (32, 64, 256)]:
        if not os.path.isdir(os.path.join(cell, cond)):
            continue
        t, uni = truth(cell, cond)
        got = {"mod": removed(os.path.join(SAS, run, f"47_spatialartifacts_mad_{b}um", f"{cond}.txt")),
               "mad": removed(os.path.join(cell, "results", mad, f"{cond}.txt")),
               "sa": removed(os.path.join(cell, "results", sa, f"{cond}.txt"))}
        for k, g in got.items():
            if g is None:
                continue
            g &= uni
            x = dict(section=sec, chem=chem, depth=dep, cond=cond, arm=k)
            if cond == "clean":
                x.update(fpr=len(g) / len(uni), band=len(g & band[0]) / len(band[0]) if band else np.nan)
            else:
                x.update(recall=len(g & t) / len(t))
            M.append(x)
        if dep == 100 and cond.endswith("_64um"):
            g45, g26, g19 = (removed(os.path.join(cell, "results", a, f"{cond}.txt"))
                             for a in (WIDE[64], AP[64], "19_spatialartifacts_8um"))
            if None not in (g45, g26, g19):
                U.append(dict(section=sec, chem=chem, cond=cond,
                              wide_sa=len((g45 | g19) & t) / len(t), ap_sa=len((g26 | g19) & t) / len(t)))
        if dep == 100 and cond in ("dropout01edge_64um", "mitoedge_64um"):
            for a in (AP[64], WIDE[64], "19_spatialartifacts_8um", "20_spatialartifacts_16um", "02_mad_8um"):
                g = removed(os.path.join(cell, "results", a, f"{cond}.txt"))
                if g is not None and chem == "probe":
                    rows.append(dict(part="D", chem=chem, depth=dep, cond=cond, a=a, section=sec,
                                     barcodes=len(g & uni & t) / len(t)))
M, U = pd.DataFrame(M), pd.DataFrame(U)

print("B. step 3 at the edge, 64 um LQRs, 100 UMIs per bin: 768 neighbours + SA 8 um vs A&P 32 um + SA 8 um")
for chem in ("probe", "3prime"):
    for cond in ("dropout01edge_64um", "mitoedge_64um", "dropout01_64um", "mito_64um"):
        u = U[(U.chem == chem) & (U.cond == cond)].set_index("section")
        if len(u):
            r = test(u.wide_sa, u.ap_sa); say(f"{chem} {cond}", r)
            rows.append(dict(part="B", chem=chem, depth=100, cond=cond, arm_a="45+19", arm_b="26+19", **r))

print("C. modified region test against the MAD threshold and SpatialArtifacts on the same bins")
for dep in (100, 12):
    for chem in ("probe", "3prime"):
        print(f" {chem}, {dep} UMIs per bin, {SAME[dep][0]} um bins")
        m = M[(M.depth == dep) & (M.chem == chem)]
        for cond in ["clean"] + [f"{k}_{dm}um" for k in KINDS for dm in (32, 64, 256)]:
            c = m[m.cond == cond]
            for val in (("fpr", "band") if cond == "clean" else ("recall",)):
                v = {k: c[c.arm == k].set_index("section")[val] for k in ("mod", "mad", "sa")}
                for other in ("mad", "sa"):
                    r = test(v["mod"], v[other]); say(f"{cond} {val}: modified vs {other}", r)
                    rows.append(dict(part="C", chem=chem, depth=dep, cond=f"{cond}:{val}", arm_a="mod", arm_b=other, **r))

print("D. changed Table 1 cells, 64 um edge rows, 100 UMIs per bin, probe-based")
d = pd.DataFrame([r for r in rows if r["part"] == "D"])
for (cond, a), g in d.groupby(["cond", "a"]):
    ref = rec(a, cond.rsplit("_", 1)[0], 64, 100, "probe").reindex(g.section.values)
    diff = np.nanmax(np.abs(g.barcodes.values - ref.values))
    print(f"  {a:26s} {cond:20s} n={len(g):2d} barcodes {100 * g.barcodes.median():5.1f}  "
          f"scores.csv {100 * ref.median():5.1f}  largest per-section difference {100 * diff:.4f} points")

out = os.path.join(GRID, "figures", "edge_checks.csv")
pd.DataFrame([r for r in rows if r["part"] != "D"]).to_csv(out, index=False)
M.to_csv(os.path.join(GRID, "figures", "edge_checks_modified_region.csv"), index=False)
print("wrote", out)
