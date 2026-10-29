import argparse, glob, os, re
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import panel as P

TOOLS = {
    "scuttle MAD":      ("02_mad_8um", "10_mad_64um"),
    "SpotSweeper":      (P.SPOTSWEEPER[8], P.SPOTSWEEPER[64]),
    "SpatialArtifacts": ("19_spatialartifacts_8um", "20_spatialartifacts_16um"),
}
SMALL_STEP = {"SpatialArtifacts": "16 µm"}
C_ALONE, C_DECONV = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d8d7d2"
DIAMS = (32.0, 64.0, 256.0)
DEPTHS = (100, 50, 25, 12)

def style():
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": INK2, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "figure.facecolor": "white", "axes.facecolor": "white"})

def load_grid(root):
    scores, markers = {d: {} for d in DEPTHS}, {d: {} for d in DEPTHS}
    for ds, dep, cell in P.cells(root, DEPTHS):
        scores[dep][ds] = P.scores(cell)
        mk = os.path.join(cell, "scored", "marker.csv")
        if os.path.exists(mk):
            markers[dep][ds] = pd.read_csv(mk)
    return scores, markers

def paired(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    k = np.isfinite(a) & np.isfinite(b)
    a, b = a[k], b[k]
    d = b - a
    d = d[~np.isclose(d, 0)]
    if len(a) < 3 or not len(d):
        return np.nan
    return stats.binomtest(int((d > 0).sum()), len(d), 0.5).pvalue

def holm(ps):
    ps = np.asarray(ps, float)
    out = np.full(len(ps), np.nan)
    idx = np.where(np.isfinite(ps))[0]
    idx = idx[np.argsort(ps[idx])]
    run = 0.0
    for r, i in enumerate(idx):
        run = max(run, min(1.0, (len(idx) - r) * ps[i]))
        out[i] = run
    return out

def stars(p):
    if not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""

def recall_pct(df, arm, kind, diam):
    q = (df.arm == arm) & (df["kind"] == kind) & (df.diameter_um == diam)
    v = df.loc[q, "recall"]
    return v.iloc[0] * 100 if len(v) else np.nan

def fp_pct(df, arm, kind, diam):
    q = (df.arm == arm) & (df["kind"] == kind) & (df.diameter_um == diam)
    r = df.loc[q]
    if not len(r):
        return np.nan
    den = r.n_total.iloc[0] - r.n_truth.iloc[0]
    return 100.0 * r.fp.iloc[0] / den if den > 0 else np.nan

def marker_ratio(mdf, arm):
    if mdf is None:
        return np.nan
    r = mdf.loc[mdf.arm == arm, "ratio_matched"]
    return r.iloc[0] if len(r) else np.nan

KEY_NATIVE, KEY_AD = "native 8 µm bins", "A&D, 64 µm bins (16 µm for SpatialArtifacts)"

def grid_figure(scores, markers, metric, title, ylab, out, annotate_marker=False,
                ylim=None, key=True):
    names = list(TOOLS)
    x = np.arange(len(names)); w = 0.34
    off_lo, off_hi = -w / 2 - 0.015, w / 2 + 0.015
    fig, axs = plt.subplots(len(DIAMS), len(DEPTHS),
                            figsize=(7.2, 5.25 if annotate_marker else 4.9), sharey="row")
    rows_out = []
    fp_top = 0.0
    for i, diam in enumerate(DIAMS):
        for j, dep in enumerate(DEPTHS):
            ax = axs[i, j]
            ds = sorted(scores.get(dep, {}))
            lo = np.array([[metric(scores[dep][d], TOOLS[t][0], KIND, diam)
                            for t in names] for d in ds], float) if ds else np.empty((0, 3))
            hi = np.array([[metric(scores[dep][d], TOOLS[t][1], KIND, diam)
                            for t in names] for d in ds], float) if ds else np.empty((0, 3))
            tops = np.zeros(len(names))
            for arr, off, col in ((lo, off_lo, C_ALONE), (hi, off_hi, C_DECONV)):
                if not arr.size:
                    continue
                m = np.nanmedian(arr, axis=0)
                q1, q3 = np.nanpercentile(arr, 25, axis=0), np.nanpercentile(arr, 75, axis=0)
                ax.bar(x + off, m, w, color=col, edgecolor="none", zorder=3)
                ax.errorbar(x + off, m, yerr=np.vstack([m - q1, q3 - m]), fmt="none", ecolor=INK2,
                            elinewidth=0.6, capsize=1.5, capthick=0.6, zorder=5)
                jit = np.linspace(-w * 0.26, w * 0.26, max(arr.shape[0], 1))
                for k in range(arr.shape[0]):
                    ax.scatter(x + off + jit[k], arr[k], s=2.4, color="#ffffff",
                               edgecolors=INK2, linewidths=0.3, zorder=6, clip_on=False)
                for xi in range(len(names)):
                    c = [q3[xi] if np.isfinite(q3[xi]) else 0.0]
                    if np.isfinite(arr[:, xi]).any():
                        c.append(np.nanmax(arr[:, xi]))
                    tops[xi] = max(tops[xi], np.nanmax(c))
                fp_top = max(fp_top, np.nanmax(tops))
            if lo.size and hi.size:
                adj = holm([paired(lo[:, xi], hi[:, xi]) for xi in range(len(names))])
                for xi, t in enumerate(names):
                    rows_out.append(dict(figure=title, lqr_um=diam, depth=dep,
                                         tool=t, n=lo.shape[0],
                                         median_alone=np.nanmedian(lo[:, xi]),
                                         median_deconv=np.nanmedian(hi[:, xi]),
                                         p_sign_holm=adj[xi]))
                    s_ = stars(adj[xi])
                    if s_:
                        ax.text(xi, tops[xi] * 1.04, s_, ha="center", va="bottom",
                                fontsize=7, color=INK, zorder=8)
            if annotate_marker:
                for xi, t in enumerate(names):
                    for slot, col, dy in ((0, C_ALONE, -0.055), (1, C_DECONV, -0.115)):
                        vals = [marker_ratio(markers.get(dep, {}).get(d), TOOLS[t][slot])
                                for d in ds]
                        vals = [v for v in vals if np.isfinite(v)]
                        ax.text(xi, dy, f"{np.median(vals):.2f}" if vals else "n.d.",
                                transform=ax.get_xaxis_transform(), ha="center",
                                va="top", fontsize=5.4, color=col)
                if i == len(DIAMS) - 1:
                    ax.tick_params(axis="x", pad=13)
            if i == 0:
                ax.set_title(f"{dep} UMIs per bin", pad=5, color=INK)
            if j == 0:
                ax.set_ylabel(f"{int(diam)} µm LQR\n{ylab}", fontsize=7,
                              fontweight="bold")
            ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0); ax.set_axisbelow(True)
            if ylim:
                ax.set_ylim(*ylim)
            ax.set_xticks(x)
            if i == len(DIAMS) - 1:
                ax.set_xticklabels(names, rotation=45, ha="right", rotation_mode="anchor",
                                   fontsize=6.2)
            else:
                ax.set_xticklabels([])
    if ylim is None:
        for ax in axs.ravel():
            ax.set_ylim(0, fp_top * (1.5 if annotate_marker else 1.25))
    if key:
        fig.legend(handles=[Patch(facecolor=C_ALONE, label=KEY_NATIVE),
                            Patch(facecolor=C_DECONV, label=KEY_AD)],
                   loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0),
                   handlelength=1.3, columnspacing=1.6)
    fig.tight_layout(w_pad=1.4, h_pad=1.0, rect=(0, 0.045 if key else 0, 1, 1))
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{out}.{ext}", bbox_inches="tight")
        print("wrote", f"{out}.{ext}")
    plt.close(fig)
    pd.DataFrame(rows_out).to_csv(f"{out}_stats.csv", index=False)

def main():
    global KIND
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="/work/benchmark_grid")
    ap.add_argument("--outdir", default="/work/benchmark_grid/figures")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    style()
    scores, markers = load_grid(a.grid)
    have = {d: len(scores[d]) for d in DEPTHS}
    print("datasets per depth:", have)
    if not any(have.values()):
        raise SystemExit("no scored grid cells found under " + a.grid)

    for ch in P.CHEMS:
        keep = lambda tbl: {d: {k: v for k, v in tbl[d].items() if P.chemistry(k) == ch}
                            for d in tbl}
        sc, mk, tag, lab = keep(scores), keep(markers), f"_{ch}", f" ({P.CHEM_LABEL[ch]})"
        last = ch == P.CHEMS[-1]
        KIND = "mito"
        grid_figure(sc, mk, recall_pct,
                    "Mitochondrial elevation: recall" + lab,
                    "Recall (%)", os.path.join(a.outdir, "fig1_1_recall_mito" + tag),
                    ylim=(0, 108), key=last)
        KIND = "dropout01"
        grid_figure(sc, mk, recall_pct,
                    "Signal loss (1% of counts retained): recall" + lab,
                    "Recall (%)", os.path.join(a.outdir, "fig1_2_recall_signal_loss" + tag),
                    ylim=(0, 108), key=last)
        KIND = "dropout01"
        grid_figure(sc, mk, fp_pct,
                    "Tissue discarded outside the injected LQRs "
                    "(upper bound on error: assumes only injected bins are damaged)" + lab,
                    "Tissue removed (%)", os.path.join(a.outdir, "fig1_3_false_positive" + tag),
                    annotate_marker=True, key=last)

if __name__ == "__main__":
    main()
