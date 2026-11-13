#!/usr/bin/env python3
import argparse, glob, os, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import panel as P
import figstyle as FS

AGG, WIDE = P.SPOTSWEEPER[64], P.WIDENED[64]
C_AGG, C_WIDE = "#eb6834", "#1baf7a"
M_AGG, M_WIDE = "o", "s"
DIAM = [32.0, 64.0, 256.0]
DEPTHS = [(100, "100 UMIs per bin"), (12, "12 UMIs per bin")]
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
    ap.add_argument("--kind", default="dropout01")
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
    fig, axs = plt.subplots(len(P.CHEMS), 2, figsize=(6.0, 5.0), sharey=True, sharex=True)
    ns = {}
    for r, ch in enumerate(P.CHEMS):
        dc = d[d.chem == ch]
        for ax, (dep, title) in zip(axs[r], DEPTHS):
            x = np.arange(len(DIAM)); w = 0.36
            for off, arm, col, mk, lab in ((-w / 2 - 0.01, AGG, C_AGG, M_AGG, "A&P, 64 µm bins"),
                                           (+w / 2 + 0.01, WIDE, C_WIDE, M_WIDE,
                                            "native 8 µm bins, 3072 neighbours")):
                both = None
                for L in DIAM:
                    for am in (AGG, WIDE):
                        secs = set(dc[(dc.depth == dep) & (dc.arm == am) & (dc["kind"] == a.kind)
                                      & (dc.diameter_um == L)].sec)
                        both = secs if both is None else both & secs
                ns[(ch, dep, lab)] = len(both)
                vals = []
                for i, L in enumerate(DIAM):
                    sub = dc[(dc.depth == dep) & (dc.arm == arm) & (dc["kind"] == a.kind)
                             & (dc.diameter_um == L) & (dc.sec.isin(both))]
                    v = sub.recall.to_numpy() * 100
                    vals.append(np.median(v) if len(v) else np.nan)
                    ax.plot(np.full(len(v), x[i] + off), v, mk, ms=3.0, mfc="white",
                            mec=col, mew=0.9, zorder=5)
                ax.bar(x + off, vals, width=w, color=col, label=lab if r == 0 else None, zorder=3)
            ax.set_xticks(x); ax.set_xticklabels([f"{int(v)} µm" for v in DIAM])
            ax.set_ylim(0, 100)
            ax.tick_params(colors=INK2)
        axs[r, 0].set_ylabel("Recall (%)", color=INK)
        axs[r, 1].text(1.04, 0.5, P.CHEM_LABEL[ch], transform=axs[r, 1].transAxes, rotation=270,
                       ha="left", va="center", fontsize=7.5, fontweight="bold", color=INK)
    fig.text(0.53, 0.075, "LQR diameter", ha="center", fontsize=7.5, color=INK)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.53, -0.01))
    fig.tight_layout(rect=(0.02, 0.1, 0.97, 1))
    for ax, letter in zip(axs[0], "ab"):
        FS.letter_axes(fig, ax, letter)
    for e in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(a.outdir, f"figS10_widened_neighbourhood.{e}"), bbox_inches="tight")
    print("wrote figS10_widened_neighbourhood; n per arm:", ns)

if __name__ == "__main__":
    main()
