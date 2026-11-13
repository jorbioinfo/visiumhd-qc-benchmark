#!/usr/bin/env python3
import argparse, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

import panel as P

COND = {"dropout01_256um": "Signal loss", "mito_256um": "Mitochondrial elevation"}
DEPTHS = ((100, "100 UMIs per bin"), (12, "12 UMIs per bin"))
TOOLS = {
    "MAD": {
        "native": ("02_mad_8um", "native 8 µm bins"),
        "ap": ("10_mad_64um", "A&P, 64 µm bins"),
        "third": ("35_mad_64um_excludepartial",
                  "A&P, 64 µm bins, partly covered bins left out")},
    "miQC": {
        "native": ("30_miqc_8um", "native 8 µm bins"),
        "ap": ("31_miqc_64um", "A&P, 64 µm bins")},
    "SpotSweeper": {
        "native": (P.SPOTSWEEPER[8], "native 8 µm bins"),
        "ap": (P.SPOTSWEEPER[64], "A&P, 64 µm bins"),
        "third": (P.WIDENED[64], "native 8 µm bins, 3072 neighbours")},
    "SpatialArtifacts": {
        "native": (P.SPATIALARTIFACTS[8], "native 8 µm bins"),
        "ap": (P.SPATIALARTIFACTS[16], "A&P, 16 µm bins")},
}

MS, DOT_ALPHA = 3.6, 0.45
CLR = {"native": "#2a78d6", "ap": "#eb6834", "third": "#1baf7a"}
ARROW_MIN_N, ARROW_MIN_RECALL, ARROW_MIN_FPR_RATIO = 5, 20.0, 4.0
INK, INK2, MUTED, BASE = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
XFLOOR, XCEIL = 0.08, 160.0
NO_CONTRAST = {("mito_256um", P.PREFIX + s) for s in P.MITO_NO_CONTRAST}
NO_RESULT = {"30_miqc_8um": "native 8 µm: mixture model does not\nconverge at this depth "
                            "(no section)"}

def style():
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
        "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8,
        "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "legend.fontsize": 6.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": BASE, "axes.linewidth": 0.6,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })

def arrow(ax, mu0, mu1, color):
    if abs(mu1[1] - mu0[1]) < ARROW_MIN_RECALL and abs(mu1[0] - mu0[0]) < np.log10(ARROW_MIN_FPR_RATIO):
        return False
    ax.annotate("", xy=(10 ** mu1[0], mu1[1]), xytext=(10 ** mu0[0], mu0[1]),
                arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.26",
                                mutation_scale=10, color=color, lw=1.6, shrinkA=0, shrinkB=0,
                                path_effects=[pe.withStroke(linewidth=3.4, foreground="white")]),
                zorder=8)
    return True

def cloud(df, arm, cond, depth, for_stats=False):
    g = df[(df.cond == cond) & (df.depth == depth) & (df.arm == arm)]
    if for_stats:
        g = g[[(cond, sec) not in NO_CONTRAST for sec in g.section]]
    if not len(g):
        return np.empty((0, 2))
    x = np.clip(g.cost.to_numpy() * 100, XFLOOR, XCEIL)
    return np.column_stack([np.log10(x), g.recall.to_numpy() * 100])

def panel(ax, df, arms, cond, depth, dlabel, first_col, first_row):
    ns = {}
    base = cloud(df, arms["native"][0], cond, depth, for_stats=True)
    for key, (arm, _lab) in arms.items():
        pts = cloud(df, arm, cond, depth)
        stat = cloud(df, arm, cond, depth, for_stats=True)
        ns[key] = len(stat)
        if not len(pts):
            continue
        if key != "native" and len(base) >= ARROW_MIN_N and len(stat) >= ARROW_MIN_N:
            arrow(ax, base.mean(0), stat.mean(0), CLR[key])
        ax.plot(10 ** pts[:, 0], pts[:, 1], "o", ms=MS, ls="none",
                zorder=5 + (key != "native"), mfc=CLR[key], mec="none", alpha=DOT_ALPHA)
    ax.axvline(50, color=BASE, lw=0.7, zorder=1)
    ax.axhline(50, color=BASE, lw=0.7, zorder=1)
    ax.axvline(100, color=BASE, lw=0.7, ls=":", zorder=1)
    for recall in (0, 100):
        ax.axhline(recall, color=BASE, lw=0.7, ls=":", zorder=1)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 1, 10, 100]); ax.set_xticklabels(["0.1", "1", "10", "100"])
    ax.set_yticks([0, 50, 100]); ax.tick_params(colors=INK2)
    if first_row:
        ax.set_title(dlabel, pad=4, color=INK)
    if not first_col:
        ax.tick_params(axis="y", labelleft=False)
    return ns

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ci", default="/work/benchmark_grid/figures/region_ci.csv")
    ap.add_argument("--outdir", default="/work/benchmark_grid/figures")
    a = ap.parse_args()
    style()
    df = pd.read_csv(a.ci)
    figures = []
    for tool, arms in TOOLS.items():
        cols = [(dep, dlabel, ch) for dep, dlabel in DEPTHS for ch in P.CHEMS]
        fig, axs = plt.subplots(2, len(cols), figsize=(7.2, 3.45), sharex=True, sharey=True)
        ns = {}
        for i, (cond, ctitle) in enumerate(COND.items()):
            for j, (dep, dlabel, ch) in enumerate(cols):
                sub = df[df.chem == ch]
                for k, v in panel(axs[i, j], sub, arms, cond, dep,
                                  f"{dlabel}\n{P.CHEM_LABEL[ch]}", j == 0, i == 0).items():
                    ns.setdefault((k, ch, dep), set()).add(v)
            axs[i, 0].set_ylabel(f"{ctitle}\nRecall (%)", fontsize=7.5, color=INK,
                                 fontweight="bold")
        fig.supxlabel("False-positive rate (%, log scale)", y=0.075, fontsize=7.5,
                      color=INK)
        hs, ls = [], []
        for key, (_arm, lab) in arms.items():
            hs.append(Line2D([], [], marker="o", ls="none", ms=MS + 1.4, mfc=CLR[key],
                             mec="none", alpha=0.75))
            ls.append(lab)
        lg = fig.legend(hs, ls, loc="lower center", ncol=len(hs), frameon=False,
                        handlelength=1.0, bbox_to_anchor=(0.5, 0.0))
        nfmt = lambda s: (f"{min(s)}" if len(s) == 1 else f"{min(s)}-{max(s)}") if s else "0"
        cov = "; ".join(P.CHEM_LABEL[ch] + " " + ", ".join(
                            f"{lab.split(',')[0].split(':')[0]} n={nfmt(ns.get((k, ch, DEPTHS[0][0]), set()))} high, "
                            f"{nfmt(ns.get((k, ch, DEPTHS[-1][0]), set()))} low"
                            for k, (_a, lab) in arms.items()) for ch in P.CHEMS)
        print(f"{tool} sections: {cov}")
        slug = tool.lower().replace(" ", "")
        figures.append((slug, fig, ns))

    x_low, x_high = zip(*(fig.axes[0].get_xlim() for _, fig, _ in figures))
    xlim = (min(XFLOOR, *x_low), max(XCEIL, *x_high))
    ylim = (-5, 105)
    for slug, fig, ns in figures:
        fig.axes[0].set(xlim=xlim, ylim=ylim)
        for ax in fig.axes:
            left, bottom, right, top = ax.dataLim.extents
            assert xlim[0] < left <= right < xlim[1] and ylim[0] < bottom <= top < ylim[1], slug
        fig.tight_layout(rect=(0, 0.1, 1, 1), h_pad=1.0, w_pad=0.7)
        for ext in ("pdf", "png", "svg"):
            fig.savefig(os.path.join(a.outdir, f"recall_fpr_{slug}.{ext}"), bbox_inches="tight")
        print("wrote", f"recall_fpr_{slug}", "n per panel:", ns)
        plt.close(fig)

if __name__ == "__main__":
    main()
