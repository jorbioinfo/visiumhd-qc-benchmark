#!/usr/bin/env python3
import argparse, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import panel as P
import figstyle as FS

PAIR = {32.0: 16, 64.0: 32, 256.0: 64}
NEIGHBOURS = {16: 192, 32: 768, 64: 3072}
C_AGG, C_WIDE = "#eb6834", "#1baf7a"
M_AGG, M_WIDE = "o", "s"
DIAM = [32.0, 64.0, 256.0]
DEPTHS = [(100, "100 UMIs per bin"), (12, "12 UMIs per bin")]
KINDS = [("dropout01", "Signal loss"), ("mito", "Mitochondrial\nelevation")]
INK, INK2, MUTED, BASE = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"

def load(grid):
    out = []
    for sec, dep, cell in P.cells(grid):
        f = P.scores(cell)
        f["sec"] = sec; f["depth"] = dep
        out.append(f)
    return pd.concat(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="/work/benchmark_grid")
    ap.add_argument("--outdir", default="/work/benchmark_grid/figures")
    a = ap.parse_args()
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600, "font.size": 7,
        "axes.labelsize": 7.5, "axes.titlesize": 8, "xtick.labelsize": 6.8,
        "ytick.labelsize": 6.8, "legend.fontsize": 6.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": BASE, "axes.linewidth": 0.6,
        "figure.facecolor": "white", "axes.facecolor": "white"})
    d = load(a.grid)
    d["chem"] = d.sec.map(P.chemistry)
    rows = [(ch, kind, ktitle) for ch in P.CHEMS for kind, ktitle in KINDS]
    fig, axs = plt.subplots(len(rows), len(DEPTHS), figsize=(6.0, 7.6), sharey=True, sharex=True)
    ns, stats = {}, []
    x = np.arange(len(DIAM)); w = 0.36
    for r, (ch, kind, ktitle) in enumerate(rows):
        dc = d[(d.chem == ch) & (d["kind"] == kind)]
        for ax, (dep, title) in zip(axs[r], DEPTHS):
            dd = dc[dc.depth == dep]
            for i, L in enumerate(DIAM):
                agg, wide = P.SPOTSWEEPER[PAIR[L]], P.WIDENED[PAIR[L]]
                sub = dd[dd.diameter_um == L]
                both = set(sub[sub.arm == agg].sec) & set(sub[sub.arm == wide].sec)
                ns[(ch, kind, dep, int(L))] = len(both)
                for off, arm, col, mk, lab in (
                        (-w / 2 - 0.01, agg, C_AGG, M_AGG, "A&P, 16, 32 and 64 µm bins"),
                        (+w / 2 + 0.01, wide, C_WIDE, M_WIDE,
                         "native 8 µm bins, 192, 768 and 3072 neighbours")):
                    v = sub[(sub.arm == arm) & sub.sec.isin(both)].recall.to_numpy() * 100
                    med = np.median(v) if len(v) else np.nan
                    ax.plot(np.full(len(v), x[i] + off), v, mk, ms=2.8, mfc="white",
                            mec=col, mew=0.9, zorder=5)
                    ax.bar(x[i] + off, med, width=w, color=col, zorder=3,
                           label=lab if (r == 0 and i == 0 and dep == DEPTHS[0][0]) else None)
                    stats.append(dict(chem=ch, kind=kind, depth=dep, lqr_um=int(L), arm=arm,
                                      n=len(v), median_recall=med))
            ax.set_xticks(x)
            ax.set_xticklabels([f"{int(L)} µm\n{PAIR[L]} µm · {NEIGHBOURS[PAIR[L]]}" for L in DIAM])
            ax.set_ylim(0, 100)
            ax.tick_params(colors=INK2)
            ax.yaxis.grid(True, color="#d8d7d2", lw=0.5, zorder=0); ax.set_axisbelow(True)
            if r == 0:
                ax.set_title(title, pad=4, color=INK)
        axs[r, 0].set_ylabel("Recall (%)", color=INK)
        axs[r, -1].text(1.04, 0.5, f"{P.CHEM_LABEL[ch]}\n{ktitle}", transform=axs[r, -1].transAxes,
                        rotation=270, ha="left", va="center", fontsize=7.5, fontweight="bold",
                        color=INK, linespacing=1.1)
    fig.supxlabel("LQR diameter, over the A&P bin size · neighbours of the pair", fontsize=7.5,
                  color=INK, y=0.045)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0.02, 0.07, 0.96, 1))
    for ax, letter in zip(axs[0], "ab"):
        FS.letter_axes(fig, ax, letter)
    for e in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(a.outdir, f"figS10_widened_neighbourhood.{e}"), bbox_inches="tight")
    pd.DataFrame(stats).to_csv(os.path.join(a.outdir, "figS10_widened_neighbourhood_stats.csv"), index=False)
    print("wrote figS10_widened_neighbourhood; sections per pair:", ns)

if __name__ == "__main__":
    main()
