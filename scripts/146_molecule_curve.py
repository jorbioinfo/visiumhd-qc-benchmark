#!/usr/bin/env python3
import argparse, math, os
from decimal import Decimal, ROUND_HALF_UP
import numpy as np, pandas as pd

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import panel as P

GRID, MOL = os.path.join(WK, "benchmark_grid"), os.path.join(WK, "benchmark_molcurve")
NEW = (200, 90, 80, 70, 60, 22, 20, 18, 15, 6, 5, 4, 3)
BINS = (8, 16, 32, 64)
TOL = 1.25
COL = {8: "#2a78d6", 16: "#1baf7a", 32: "#eb6834", 64: "#7a5cc7"}
INK, INK2, GRIDC = "#0b0b0b", "#52514e", "#d8d7d2"
NAME = {P.PROBE: "probe-based", P.POLYA: "3′ polyA"}

def sign_p(up, dn):
    n = up + dn
    return float("nan") if n == 0 else min(1.0, 2 * sum(math.comb(n, k) for k in range(min(up, dn) + 1)) / 2 ** n)

def fp(p):
    if p != p:
        return "n.a."
    if p < 0.001:
        return "< 0.001"
    if p >= 1:
        return "1.0"
    d = Decimal(p)
    return str(d.quantize(Decimal(1).scaleb(d.adjusted() - 1), rounding=ROUND_HALF_UP))

def load():
    arm = {v: k for k, v in P.MAD.items() if k in BINS}
    missing = set(NEW) - {depth for _, depth, _ in P.cells(MOL, NEW)}
    if missing:
        raise SystemExit(f"Missing scored MAD runs in {MOL} at depths {sorted(missing)}. See the README for Figure S9 inputs.")
    rows = []
    for root, deps in ((GRID, P.DEPTHS), (MOL, NEW)):
        for sec, dep, cell in P.cells(root, deps):
            f = pd.read_csv(os.path.join(cell, "scored", "scores.csv"))
            f = f[(f.condition == "clean") & f.arm.isin(arm)]
            for _, r in f.iterrows():
                rows.append(dict(section=sec, chem=P.chemistry(sec), depth=dep, bin=arm[r.arm], fpr=r.discard_frac))
    D = pd.DataFrame(rows).groupby(["chem", "section", "bin", "depth"], as_index=False).fpr.mean()
    D["mol"] = D.depth * (D["bin"] / 8) ** 2
    ref = D[(D["bin"] == 8) & (D.depth == 100)].set_index("section").fpr
    own = D[D.depth == 100].set_index(["section", "bin"]).fpr
    D = D[D.section.isin(ref.index)].copy()
    D["ratio8"] = D.fpr / D.section.map(ref)
    D["ratio"] = D.fpr / [own.get((s, b), np.nan) for s, b in zip(D.section, D["bin"])]
    return D, ref

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=os.path.join(GRID, "figures"))
    args = parser.parse_args()
    D, ref = load()
    lines = [f"MAD threshold, uninjected condition. Ratio = false-positive rate / the section's rate on the same "
             f"bin size at 100 UMIs per bin; ratio to 8 µm = / its rate on 8 µm bins at 100. Sign tests against the "
             f"same bin size at 100 and against 8 µm bins at 100, ties left out.", ""]
    for ch in P.CHEMS:
        d = D[D.chem == ch]
        lines.append(f"## {NAME[ch]} (sections with 100 UMIs per bin: {d.section.nunique()})")
        lines.append(f"{'molecules':>9} {'bin':>5} {'UMIs/bin':>8} {'n':>3} {'median FPR %':>12} {'ratio':>6} "
                     f"{'above 1.25x':>11}  {'sign test, same bin at 100':<32} {'ratio to 8 µm':>13}  sign test, 8 µm at 100")
        for (b, dep), g in sorted(d.groupby(["bin", "depth"]), key=lambda x: (x[0][0] * x[0][0] * x[0][1], x[0][0])):
            if dep == 100:
                continue
            r = g.dropna(subset=["ratio"])
            u1, d1 = int((r.ratio > 1).sum()), int((r.ratio < 1).sum())
            u8, d8 = int((g.ratio8 > 1).sum()), int((g.ratio8 < 1).sum())
            lines.append(f"{dep * (b / 8) ** 2:>9g} {b:>3} µm {dep:>8g} {len(g):>3} {100 * g.fpr.median():>12.1f} "
                         f"{r.ratio.median():>6.2f} {int((r.ratio > TOL).sum()):>5}/{len(r):<5}  "
                         f"{f'higher {u1}, lower {d1}; p {fp(sign_p(u1, d1))}':<32} {g.ratio8.median():>13.2f}  "
                         f"higher {u8}, lower {d8}; p {fp(sign_p(u8, d8))}")
        lines.append("")
    lines.append(f"## Per-section threshold: fewest molecules per bin from which every bin size and depth stays within "
                 f"{TOL:g} times the section's rate on the same bin size at 100 UMIs per bin")
    th = {}
    for (ch, sec), g in D.groupby(["chem", "section"]):
        g = g[(g.depth != 100) & g.ratio.notna()]
        ms = sorted(g.mol.unique())
        t = next((m for m in ms if (g[g.mol >= m].ratio <= TOL).all()), float("nan"))
        th.setdefault(ch, []).append(t)
    for ch in P.CHEMS:
        v = np.array(th.get(ch, []), float)
        if len(v):
            q = np.nanpercentile(v, [25, 50, 75])
            ok = v[~np.isnan(v)]
            lines.append(f"{NAME[ch]}: n = {len(v)}, median {q[1]:g}, interquartile range {q[0]:g}–{q[2]:g}, "
                         f"max {np.nanmax(v):g}, no threshold {int(np.isnan(v).sum())}; "
                         + ", ".join(f"{m:g}: {int((ok == m).sum())}" for m in sorted(set(ok))))
    lines.append("")
    med = D[D.depth != 100].groupby(["chem", "bin", "depth"]).agg(mol=("mol", "first"),
                                                                                         ratio=("ratio", "median"))
    ms = sorted(med.mol.unique())
    pooled = next(m for m in ms if (med[med.mol >= m].ratio <= TOL).all())
    lines.append(f"Pooled: from {pooled:g} molecules per bin, every bin size and depth has a median ratio within "
                 f"{TOL:g} in both chemistries.")
    os.makedirs(args.outdir, exist_ok=True)
    open(os.path.join(args.outdir, "molecule_curve.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    draw(D, args.outdir)

def draw(D, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600, "font.size": 7,
                         "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 7,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2,
                         "axes.linewidth": 0.6, "figure.facecolor": "white", "axes.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=True)
    for ax, ch, letter in zip(axes, P.CHEMS, "ab"):
        d = D[D.chem == ch]
        for b in BINS:
            g = d[d["bin"] == b].groupby("depth").agg(mol=("mol", "first"), med=("ratio", "median"),
                                                     lo=("ratio", lambda x: np.percentile(x, 25)),
                                                     hi=("ratio", lambda x: np.percentile(x, 75))).sort_values("mol")
            ax.errorbar(g.mol, g.med, yerr=[g.med - g.lo, g.hi - g.med], color=COL[b], marker="o", ms=3, lw=1,
                        elinewidth=0.6, capsize=0, label=f"{b} µm bins")
        ax.axhline(1, color=GRIDC, lw=0.8, zorder=0)
        ax.axhline(TOL, color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=0)
        ax.axvline(100, color=GRIDC, lw=0.8, zorder=0)
        ax.set_xscale("log")
        ticks = [3, 6, 12, 25, 50, 100, 200, 800, 3200, 12800]
        ax.set_xticks(ticks); ax.set_xticklabels([f"{x:,}" for x in ticks]); ax.minorticks_off()
        ax.set_xlabel("Molecules per bin (UMIs per 8 µm bin × 8 µm bins summed)")
        ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", color=INK,
                family="DejaVu Sans", va="bottom")
    axes[0].set_ylabel("False-positive rate relative to\nthe same bin size at 100 UMIs per bin")
    axes[1].legend(frameon=False, loc="upper right")
    fig.tight_layout()
    base = os.path.join(outdir, "figS9_molecule_curve")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{base}.{ext}")
    print("wrote", base + ".{png,pdf,svg}")

if __name__ == "__main__":
    main()
