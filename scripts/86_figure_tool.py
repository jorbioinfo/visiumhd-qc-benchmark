import argparse, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import panel as P

COND = {"dropout01_256um": "Signal loss", "mito_256um": "Mitochondrial elevation"}
DEPTHS = ((100, "100 UMIs per bin"), (12, "12 UMIs per bin"))
TOOLS = {
    "MAD": {
        "o": ("02_mad_8um", "native 8 µm bins"),
        "*": ("10_mad_64um", "A&D, 64 µm bins"),
        "^": ("35_mad_64um_excludepartial",
              "A&D, 64 µm bins, partly covered bins left out")},
    "miQC": {
        "o": ("30_miqc_8um", "native 8 µm bins"),
        "*": ("31_miqc_64um", "A&D, 64 µm bins")},
    "SpotSweeper": {
        "o": (P.SPOTSWEEPER[8], "native 8 µm bins"),
        "*": (P.SPOTSWEEPER[64], "A&D, 64 µm bins"),
        "s": (P.WIDENED[64], "native 8 µm bins, 3072 neighbours")},
}

MARK = {"o": ("o", 5.4, False), "*": ("o", 6.0, True),
        "s": ("s", 5.6, True), "^": ("^", 6.4, True)}
CLR = {"o": "#2a78d6", "*": "#eb6834", "s": "#1baf7a", "^": "#1baf7a"}
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

ELLIPSE_MIN_N = 5

def mean_ellipse(ax, pts, color):
    n = len(pts)
    if n < ELLIPSE_MIN_N:
        return
    mu = pts.mean(0)
    S = np.cov(pts.T, ddof=1)
    if not np.all(np.isfinite(S)) or np.linalg.det(S) <= 0:
        return
    scale = 2 * (n - 1) / (n * (n - 2)) * stats.f.ppf(0.95, 2, n - 2)
    w, V = np.linalg.eigh(S * scale)
    w = np.clip(w, 0, None)
    t = np.linspace(0, 2 * np.pi, 128)
    circ = np.stack([np.cos(t), np.sin(t)])
    xy = (V @ (np.sqrt(w)[:, None] * circ)).T + mu
    ax.fill(10 ** xy[:, 0], xy[:, 1], color=color, alpha=0.14, lw=0, zorder=2)
    ax.plot(10 ** xy[:, 0], xy[:, 1], color=color, lw=0.7, alpha=.5, zorder=2)

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
    base = cloud(df, arms["o"][0], cond, depth, for_stats=True)
    for key, (arm, _lab) in arms.items():
        pts = cloud(df, arm, cond, depth)
        stat = cloud(df, arm, cond, depth, for_stats=True)
        ns[key] = len(stat)
        if not len(pts):
            continue
        mean_ellipse(ax, stat, CLR[key])
        if key != "o" and len(base) and len(stat):
            mu0, mu1 = base.mean(0), stat.mean(0)
            ax.annotate("", xy=(10 ** mu1[0], mu1[1]), xytext=(10 ** mu0[0], mu0[1]),
                        arrowprops=dict(arrowstyle="-|>", color=CLR[key], lw=1.4,
                                        shrinkA=6, shrinkB=7, alpha=.85), zorder=4)
        m, ms, filled = MARK[key]
        c = CLR[key]
        ax.plot(10 ** pts[:, 0], pts[:, 1], m, ms=ms, ls="none", zorder=5 + (key != "o"),
                mfc=c if filled else "white",
                mec="white" if filled else c, mew=0.9 if filled else 1.3, alpha=.9)
    ax.axvline(50, color=BASE, lw=0.7, zorder=1)
    ax.axhline(50, color=BASE, lw=0.7, zorder=1)
    ax.set_xscale("log"); ax.set_xlim(XFLOOR, XCEIL); ax.set_ylim(-4, 106)
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
        fig.supxlabel("Undamaged tissue removed (%, log scale)", y=0.075, fontsize=7.5,
                      color=INK)
        hs, ls = [], []
        for key, (_arm, lab) in arms.items():
            m, ms, filled = MARK[key]
            c = CLR[key]
            hs.append(Line2D([], [], marker=m, ls="none", ms=ms,
                             mfc=c if filled else "white",
                             mec="white" if filled else c,
                             mew=0.9 if filled else 1.3))
            ls.append(lab)
        lg = fig.legend(hs, ls, loc="lower center", ncol=len(hs), frameon=False,
                        handlelength=1.0, bbox_to_anchor=(0.5, 0.0))
        nfmt = lambda s: (f"{min(s)}" if len(s) == 1 else f"{min(s)}-{max(s)}") if s else "0"
        cov = "; ".join(P.CHEM_LABEL[ch] + " " + ", ".join(
                            f"{lab.split(',')[0].split(':')[0]} n={nfmt(ns.get((k, ch, DEPTHS[0][0]), set()))} high, "
                            f"{nfmt(ns.get((k, ch, DEPTHS[-1][0]), set()))} low"
                            for k, (_a, lab) in arms.items()) for ch in P.CHEMS)
        print(f"{tool} sections: {cov}")
        fig.tight_layout(rect=(0, 0.1, 1, 1), h_pad=1.0, w_pad=0.7)
        slug = tool.lower().replace(" ", "")
        for ext in ("pdf", "png", "svg"):
            fig.savefig(os.path.join(a.outdir, f"fig1_{slug}.{ext}"), bbox_inches="tight")
        print("wrote", f"fig1_{slug}", "n per panel:", ns)
        plt.close(fig)

if __name__ == "__main__":
    main()
