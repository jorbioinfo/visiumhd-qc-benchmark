import argparse, glob, os, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import panel as P
import figstyle as FS

BINS = (8, 16, 32, 64)
ARMS = {"scuttle MAD": P.MAD, "SpotSweeper": P.SPOTSWEEPER}
KINDS = [("mito", "Mitochondrial elevation"), ("dropout01", "Signal loss")]
DIAMS = (32.0, 64.0, 256.0)
DIAM_COL = {32.0: "#2a78d6", 64.0: "#1baf7a", 256.0: "#eb6834"}
DIAM_MK = {32.0: "o", 64.0: "s", 256.0: "^"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d8d7d2"

def style():
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
        "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": INK2, "axes.linewidth": 0.6,
        "figure.facecolor": "white", "axes.facecolor": "white"})

def load(root, depth):
    return {sec: P.scores(cell)
            for sec, _d, cell in P.cells(root, (depth,))}

METRIC_LABEL = {"recall": "Recall (%)", "f1": "F1 (%)"}

def score(df, arm, kind, diam, metric):
    q = (df.arm == arm) & (df["kind"] == kind) & (df.diameter_um == diam)
    v = df.loc[q, metric]
    return v.iloc[0] * 100 if len(v) else np.nan

def block(sf, frames, depth, metric, drop=(), header=True, rows=None):
    cols = [(t, ch) for t in ARMS for ch in P.CHEMS]
    axs = sf.subplots(len(KINDS), len(cols), sharex=True, sharey="row")
    for i, (kind, ktitle) in enumerate(KINDS):
        for j, (t, ch) in enumerate(cols):
            ds = [d for d in frames
                  if P.chemistry(d) == ch and not any(k in d for k in drop)]
            ax = axs[i, j]
            for diam in DIAMS:
                per = np.array([[score(frames[d], ARMS[t][b], kind, diam, metric)
                                 for b in BINS] for d in ds], float)
                m = np.nanmedian(per, axis=0)
                q1, q3 = np.nanpercentile(per, 25, axis=0), np.nanpercentile(per, 75, axis=0)
                col = DIAM_COL[diam]
                ax.fill_between(BINS, q1, q3, color=col, alpha=.12, lw=0)
                ax.plot(BINS, m, color=col, lw=1.8, marker=DIAM_MK[diam], ms=4.5,
                        mec="white", mew=.8, zorder=5,
                        label=f"{int(diam)} µm LQR")
                if rows is not None:
                    for b, v in zip(BINS, m):
                        rows.append(dict(depth=depth, chem=ch, damage=kind, tool=t,
                                         lesion_um=diam, qc_bin_um=b, metric=metric, value=v,
                                         n=per.shape[0]))
            ax.set_xscale("log", base=2)
            ax.set_xticks(BINS)
            ax.set_xticklabels([str(b) for b in BINS])
            ax.set_ylim(0, 105)
            ax.yaxis.grid(True, color=GRID, lw=.5, zorder=0); ax.set_axisbelow(True)
            if i == 0 and header:
                ax.set_title(f"{t}\n{P.CHEM_LABEL[ch]}", pad=6, color=INK)
            if j == 0:
                ax.set_ylabel(METRIC_LABEL[metric], fontsize=7.5)
            if j == len(cols) - 1:
                ax.text(1.04, 0.5, ktitle.replace(" ", "\n", 1), transform=ax.transAxes,
                        rotation=270, ha="left", va="center", fontsize=7.5,
                        fontweight="bold", color=INK, linespacing=1.1)
    return axs

def save(fig, out):
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{out}.{ext}", bbox_inches="tight")
        print("wrote", f"{out}.{ext}")
    plt.close(fig)

def key(fig, axs):
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside upper center", ncol=3, frameon=False, handlelength=1.8,
               columnspacing=2.2, fontsize=7.2)

def figure(frames, out, depth, drop=(), metric="recall"):
    fig = plt.figure(figsize=(7.2, 4.4), layout="constrained")
    rows = []
    axs = block(fig, frames, depth, metric, drop, rows=rows)
    fig.supxlabel("bin size QC was run on (µm)", fontsize=7.5)
    key(fig, axs)
    save(fig, out)
    pd.DataFrame(rows).to_csv(f"{out}_stats.csv", index=False)

def stacked(grid, out, depths, metric, titles, height=3.4):
    W = 7.2
    fig = plt.figure(figsize=(W, height * len(depths) + 0.35), layout="constrained")
    sfs = fig.subfigures(len(depths), 1, hspace=0.04)
    for n, (sf, depth) in enumerate(zip(sfs, depths)):
        axs = block(sf, load(grid, depth), depth, metric, header=(n == 0))
        if titles:
            sf.suptitle(f"{depth} UMIs per bin", fontsize=8, color=INK)
        FS.panel_letter(sf, 0.0, 1.0, "abcd"[n], W)
        if n == 0:
            key(fig, axs)
        if n == len(depths) - 1:
            sf.supxlabel("bin size QC was run on (µm)", fontsize=7.5)
    save(fig, out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="/work/benchmark_grid")
    ap.add_argument("--outdir", default="/work/benchmark_grid/figures")
    ap.add_argument("--depth", type=int, default=100)
    ap.add_argument("--metric", default="recall", choices=("recall", "f1"))
    ap.add_argument("--stack", type=int, nargs="+",
                    help="draw these depths stacked into one figure instead")
    ap.add_argument("--name", help="file name of the stacked figure")
    ap.add_argument("--depth-titles", action="store_true")
    ap.add_argument("--block-height", type=float, default=3.4, help="inches per depth")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    style()
    if a.stack:
        stacked(a.grid, os.path.join(a.outdir, a.name), a.stack, a.metric, a.depth_titles,
                a.block_height)
        return
    frames = load(a.grid, a.depth)
    suffix = "" if a.metric == "recall" else f"_{a.metric}"
    figure(frames, os.path.join(a.outdir, f"fig5_binsize_d{a.depth}{suffix}"), a.depth,
           metric=a.metric)
    figure(frames, os.path.join(a.outdir, f"fig5_binsize_d{a.depth}{suffix}_noKidney"),
           a.depth, drop=tuple(P.MITO_NO_CONTRAST), metric=a.metric)

if __name__ == "__main__":
    main()
