#!/usr/bin/env python3
import math
import os
from decimal import Decimal, ROUND_HALF_UP

import numpy as np
import pandas as pd

import panel as P

WK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(WK, "benchmark_grid", "figures", "table3_tests.txt")
BINS = (8, 16, 32, 64)
DIAMS = (32, 64, 256)
KINDS = (("dropout01", "signal loss"), ("mito", "mito. elevation"))
TESTS = (("Neighbourhood test (SpotSweeper, 48 neighbours)", P.SPOTSWEEPER),
         ("Section-wide test (MAD threshold)", P.MAD))
SA = P.SPATIALARTIFACTS
SA_KINDS = (("dropout01", "inside the tissue"), ("dropout01edge", "at the tissue edge"))
FPR_ROWS = (("neighbourhood", P.SPOTSWEEPER), ("section-wide (MAD)", P.MAD), ("region (SA)", SA))

def sign_p(up, dn):
    n = up + dn
    if n == 0:
        return float("nan")
    return min(1.0, 2 * sum(math.comb(n, k) for k in range(min(up, dn) + 1)) / 2 ** n)

def fp(p):
    if p != p:
        return "n.a."
    if p < 0.001:
        return "< 0.001"
    if p >= 1:
        return "1.0"
    d = Decimal(p)
    return str(d.quantize(Decimal(1).scaleb(d.adjusted() - 1), rounding=ROUND_HALF_UP))

def count(pairs):
    up = sum(a > b for a, b in pairs)
    dn = sum(a < b for a, b in pairs)
    return up, dn, len(pairs) - up - dn

def load():
    rows = []
    for sec, dep, cell in P.cells(os.path.join(WK, "benchmark_grid")):
        f = P.scores(cell)
        f = f[f.arm.isin(list(P.SPOTSWEEPER.values()) + list(P.MAD.values()) + list(SA.values()))]
        f = f.assign(section=P.short(sec), depth=dep, chem=P.chemistry(sec))
        rows.append(f)
    s = pd.concat(rows, ignore_index=True)
    s["kind"] = s["kind"].astype(str)
    return s

def optimum(values):
    v = np.array(values, float)
    return None if np.all(np.isnan(v)) else BINS[int(np.nanargmax(v))]

def main():
    S = load()
    rec = S[S.kind.isin([k for k, _ in KINDS] + [k for k, _ in SA_KINDS]) & S.diameter_um.isin(DIAMS)]
    R = rec.pivot_table(index=["chem", "section", "kind", "diameter_um", "depth"], columns="arm",
                        values="recall")
    clean = S[S.condition == "clean"].pivot_table(index=["chem", "section", "depth"], columns="arm",
                                                  values="discard_frac")
    lines = []
    say = lines.append
    for ch in P.CHEMS:
        say("#" * 78)
        say(f"# {P.CHEM_LABEL[ch]} sections")
        say("#" * 78)
        secs = sorted(R.loc[ch].index.get_level_values("section").unique())
        for title, arms in TESTS:
            say("")
            say(f"{title}: best bin size (µm) per section, recall averaged over the depths")
            head = "".join(f"{lab[:4]} {d:>3}" .rjust(10) for _, lab in KINDS for d in DIAMS)
            say(f"{'section':42}{head}")
            pooled = {}
            for sec in secs:
                cells = []
                for k, _ in KINDS:
                    for d in DIAMS:
                        try:
                            x = R.loc[(ch, sec, k, d)]
                            o = optimum([x[arms[b]].mean() if arms[b] in x else np.nan for b in BINS])
                        except KeyError:
                            o = None
                        pooled[(sec, k, d)] = o
                        cells.append("-" if o is None else str(o))
                say(f"{sec:42}" + "".join(c.rjust(10) for c in cells))
            say("  Best bin size larger for a 256 µm-diameter LQR than for a 32 µm one:")
            for k, lab in KINDS:
                pr = [(pooled[(s, k, 256)], pooled[(s, k, 32)]) for s in secs
                      if pooled.get((s, k, 256)) and pooled.get((s, k, 32))]
                up, dn, eq = count(pr)
                say(f"    {lab:16} larger in {up}, smaller in {dn}, tied in {eq} of {len(pr)} sections; "
                    f"p {fp(sign_p(up, dn))}")
            if arms is P.SPOTSWEEPER:
                say("  Best bin size larger for mitochondrial elevation than for signal loss:")
                for d in (32, 64):
                    pr = [(pooled[(s, 'mito', d)], pooled[(s, 'dropout01', d)]) for s in secs
                          if pooled.get((s, 'mito', d)) and pooled.get((s, 'dropout01', d))]
                    up, dn, eq = count(pr)
                    say(f"    {d} µm LQRs      larger in {up}, smaller in {dn}, tied in {eq} of {len(pr)} "
                        f"sections; p {fp(sign_p(up, dn))}")
                say("  Best bin size larger at 12 than at 100 UMIs per bin (single depths):")
                for k, lab in KINDS:
                    for d in (32, 64, 256):
                        pr = []
                        for s in secs:
                            o = {}
                            for dep in (100, 12):
                                try:
                                    x = R.loc[(ch, s, k, d, dep)]
                                    o[dep] = optimum([x.get(arms[b], np.nan) for b in BINS])
                                except KeyError:
                                    pass
                            if o.get(100) and o.get(12):
                                pr.append((o[12], o[100]))
                        up, dn, eq = count(pr)
                        say(f"    {lab:16} {d:>3} µm: larger in {up}, smaller in {dn}, tied in {eq} of "
                            f"{len(pr)} sections; p {fp(sign_p(up, dn))}")
                say("  32 µm-diameter mitochondrial LQRs at 12 UMIs per bin (single depth):")
                best, pr = [], []
                for s in secs:
                    try:
                        x = R.loc[(ch, s, "mito", 32, 12)]
                    except KeyError:
                        continue
                    v = [x.get(arms[b], np.nan) for b in BINS]
                    o = optimum(v)
                    if o:
                        best.append(o)
                        pr.append((v[3], v[2]))
                up, dn, eq = count(pr)
                say("    best bin size: " + ", ".join(f"{b} µm in {best.count(b)}" for b in BINS)
                    + f" of {len(best)} sections")
                say(f"    recall on 64 against 32 µm bins: median {100 * np.median([q for q, _ in pr]):.1f} against "
                    f"{100 * np.median([q for _, q in pr]):.1f}; higher in {up}, lower in {dn}, tied in {eq} of "
                    f"{len(pr)} sections; p {fp(sign_p(up, dn))}")
        say("")
        say("Region test (SpatialArtifacts): recall (%) of signal-loss LQRs, median over sections")

        def sa_rec(arm, kind, d, dep):
            try:
                return R.loc[(ch, slice(None), kind, d, dep), arm].droplevel(
                    ["chem", "kind", "diameter_um", "depth"]).dropna()
            except KeyError:
                return pd.Series(dtype=float)
        for kind, place in SA_KINDS:
            say(f"  LQRs {place}")
            for b in (8, 16):
                for d in DIAMS:
                    r = {dep: sa_rec(SA[b], kind, d, dep) for dep in P.DEPTHS}
                    med = "".join(f"{100 * r[dep].median():>7.1f} (n={len(r[dep]):>2})" for dep in P.DEPTHS)
                    pr = pd.concat([r[12], r[100]], axis=1, keys=["lo", "hi"]).dropna()
                    up, dn, eq = count(list(zip(pr.lo, pr.hi)))
                    z = r[12]
                    say(f"    {b:>2} µm bins, {d:>3} µm LQRs: {med}; 12 against 100 UMIs/bin: higher {up}, lower {dn}, "
                        f"tied {eq} of {len(pr)}; p {fp(sign_p(up, dn))}; at 12, below 1% in {int((z < 0.01).sum())}, "
                        f"at least 50% in {int((z >= 0.5).sum())} of {len(z)} sections")
            say(f"  16 µm bins with A&P against 8 µm bins, LQRs {place}")
            for d in DIAMS:
                for dep in (100, 12):
                    a, b8 = sa_rec(SA[16], kind, d, dep), sa_rec(SA[8], kind, d, dep)
                    pr = pd.concat([a, b8], axis=1, keys=["a", "n"]).dropna()
                    up, dn, eq = count(list(zip(pr.a, pr.n)))
                    ma, mn = pr.a >= 0.5, pr.n >= 0.5
                    only_a, only_n = int((ma & ~mn).sum()), int((~ma & mn).sum())
                    say(f"    {d:>3} µm LQRs at {dep:>3} UMIs/bin: {100 * pr.a.median():.1f} against {100 * pr.n.median():.1f}; "
                        f"higher {up}, lower {dn}, tied {eq} of {len(pr)}; p {fp(sign_p(up, dn))}; at least 50% recalled in "
                        f"{int(ma.sum())} against {int(mn.sum())} sections, only on 16 µm in {only_a}, only on 8 µm in "
                        f"{only_n}; p {fp(sign_p(only_a, only_n))}")
        say("")
        say("False-positive rate (%) on the uninjected condition, median over sections")
        say(f"{'test':22}{'bin':>6}" + "".join(f"{str(d) + ' UMIs/bin':>14}" for d in P.DEPTHS)
            + f"{'12 vs 100, within a section':>44}")
        for short, arms in FPR_ROWS:
            for b in BINS:
                if arms.get(b) not in clean.columns:
                    continue
                c = clean.loc[ch, arms[b]].unstack("depth")
                med = "".join(f"{100 * c[d].median():>9.1f} (n={c[d].notna().sum():>2})"
                              if d in c else f"{'-':>14}" for d in P.DEPTHS)
                pr = c[[100, 12]].dropna()
                up, dn, eq = count(list(zip(pr[12], pr[100])))
                over12 = int((pr[12] > 1.25 * pr[100]).sum())
                p50 = c[[100, 50]].dropna()
                over50 = int((p50[50] > 1.25 * p50[100]).sum())
                none = "".join(f" {d}: {int((c[d] == 0).sum())}" for d in P.DEPTHS if d in c)
                say(f"{short:22}{b:>4} µm{med}   higher {up}, lower {dn}, tied {eq} of {len(pr)}; "
                    f"p {fp(sign_p(up, dn))}; above 1.25 times at 50: {over50}/{len(p50)}, at 12: "
                    f"{over12}/{len(pr)}" + (f"; sections with nothing removed,{none}" if arms is SA else ""))
        c16 = clean.loc[ch, P.MAD[16]].unstack("depth")[12]
        c8 = clean.loc[ch, P.MAD[8]].unstack("depth")[50]
        pr = pd.concat([c16, c8], axis=1, keys=["b16", "b8"]).dropna()
        up, dn, eq = count(list(zip(pr.b16, pr.b8)))
        say(f"section-wide (MAD), 16 µm bins at 12 against 8 µm bins at 50 UMIs/bin (about 50 molecules per bin): "
            f"{100 * pr.b16.median():.1f} against {100 * pr.b8.median():.1f}; higher {up}, lower {dn}, tied {eq} of "
            f"{len(pr)}; p {fp(sign_p(up, dn))}")
        ref = clean.loc[ch, P.MAD[8]].unstack("depth")[100]
        for b, d in ((16, 25), (16, 50), (32, 12)):
            cb = clean.loc[ch, P.MAD[b]].unstack("depth")[d]
            pr = pd.concat([cb, ref], axis=1, keys=["b", "ref"]).dropna()
            up, dn, eq = count(list(zip(pr.b, pr.ref)))
            say(f"section-wide (MAD), {b} µm bins at {d} against 8 µm bins at 100 UMIs/bin (about "
                f"{round(d * (b / 8) ** 2)} against 100 molecules per bin): {100 * pr.b.median():.1f} against "
                f"{100 * pr.ref.median():.1f}; higher {up}, lower {dn}, tied {eq} of {len(pr)}; p {fp(sign_p(up, dn))}")
        say("")
    open(OUT, "w").write("\n".join(lines) + "\n")
    print("wrote", OUT, len(lines), "lines")

if __name__ == "__main__":
    main()
