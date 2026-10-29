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

def discard_pct(df, arm):
    v = df.loc[(df.arm == arm) & (df.condition == "clean"), "discard_frac"]
    return v.iloc[0] * 100 if len(v) else np.nan

def retention(mdf, arm):
    if mdf is None:
        return np.nan
    r = mdf.loc[mdf.arm == arm, "ratio_matched"]
    return r.iloc[0] if len(r) else np.nan

def match_scale(mdf, arm):
    if mdf is None or "match_scale_um" not in mdf.columns:
        return np.nan
    r = mdf.loc[mdf.arm == arm, "match_scale_um"]
    return r.iloc[0] if len(r) else np.nan

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

def sign_p(d):
    d = np.asarray(d, float)
    d = d[np.isfinite(d) & ~np.isclose(d, 0)]
    return stats.binomtest(int((d > 0).sum()), len(d), 0.5).pvalue if len(d) else np.nan

def collect(getter, table, dep, names, slot):
    ds = sorted(table.get(dep, {}))
    if not ds:
        return np.empty((0, len(names))), []
    return np.array([[getter(table[dep][d], TOOLS[t][slot]) for t in names]
                     for d in ds], float), ds

def draw_bars(ax, arr, x, off, w, col, floor=0.0):
    if not arr.size:
        return np.zeros(len(x))
    m = np.nanmedian(arr, axis=0)
    q1, q3 = np.nanpercentile(arr, 25, axis=0), np.nanpercentile(arr, 75, axis=0)
    ax.bar(x + off, m - floor, w, bottom=floor, color=col, edgecolor="none", zorder=3)
    ax.errorbar(x + off, m, yerr=np.vstack([m - q1, q3 - m]), fmt="none", ecolor=INK2, elinewidth=0.6,
                capsize=1.5, capthick=0.6, zorder=5)
    jit = np.linspace(-w * 0.26, w * 0.26, max(arr.shape[0], 1))
    for k in range(arr.shape[0]):
        ax.scatter(x + off + jit[k], arr[k], s=2.4, color="#ffffff",
                   edgecolors=INK2, linewidths=0.3, zorder=6)
    return np.nan_to_num(np.nanmax(arr, axis=0))

KEY_NATIVE, KEY_AD = "native 8 µm bins", "A&D, 64 µm bins (16 µm for SpatialArtifacts)"

def figure(scores, markers, out, key=True):
    names = list(TOOLS)
    x = np.arange(len(names)); w = 0.34
    off_lo, off_hi = -w / 2 - 0.015, w / 2 + 0.015
    fig, axs = plt.subplots(2, len(DEPTHS), figsize=(7.2, 4.6), sharey="row")
    rows = []
    top_a, ret_lo, ret_hi = 0.0, 1.0, 1.0
    for j, dep in enumerate(DEPTHS):
        lo_d, ds = collect(discard_pct, scores, dep, names, 0)
        hi_d, _ = collect(discard_pct, scores, dep, names, 1)
        lo_r, _ = collect(retention, markers, dep, names, 0)
        hi_r, _ = collect(retention, markers, dep, names, 1)

        ax = axs[0, j]
        for arr, off, col in ((lo_d, off_lo, C_ALONE), (hi_d, off_hi, C_DECONV)):
            top_a = max(top_a, np.nanmax(draw_bars(ax, arr, x, off, w, col), initial=0.0))
        if lo_d.size and hi_d.size and lo_d.shape[0] >= 3:
            adj = holm([sign_p(hi_d[:, i] - lo_d[:, i]) for i in range(len(names))])
            for i in range(len(names)):
                s_ = stars(adj[i])
                if s_:
                    ax.text(i, max(np.nanmax(lo_d[:, i], initial=0),
                                   np.nanmax(hi_d[:, i], initial=0)) * 1.05, s_,
                            ha="center", va="bottom", fontsize=7, color=INK, zorder=8)
        if j == 0:
            ax.set_ylabel("Tissue removed (%)", fontweight="bold", labelpad=4)
        ax.set_title(f"{dep} UMIs per bin", pad=5, color=INK)

        ax = axs[1, j]
        ax.axhline(0.0, color=INK2, lw=0.7, zorder=4)
        for arr, off, col in ((lo_r, off_lo, C_ALONE), (hi_r, off_hi, C_DECONV)):
            with np.errstate(divide="ignore", invalid="ignore"):
                lg_r = np.log2(np.where(np.isfinite(arr) & (arr > 0), arr, np.nan))
            draw_bars(ax, lg_r, x, off, w, col, floor=0.0)
            if lg_r.size and np.isfinite(lg_r).any():
                q1 = np.nanpercentile(lg_r, 25, axis=0); q3 = np.nanpercentile(lg_r, 75, axis=0)
                ret_lo = min(ret_lo, np.nanmin(lg_r), np.nanmin(q1))
                ret_hi = max(ret_hi, np.nanmax(lg_r), np.nanmax(q3))
        for arr, slot in ((lo_r, 0), (hi_r, 1)):
            if arr.shape[0] < 3:
                continue
            for i, t in enumerate(names):
                v = arr[:, i][np.isfinite(arr[:, i])]
                p = sign_p(v - 1.0) if len(v) >= 3 else np.nan
                rows.append(dict(depth=dep, tool=t,
                                 arm=TOOLS[t][slot], n_datasets=len(v),
                                 discard_pct=np.nanmedian((lo_d, hi_d)[slot][:, i]),
                                 retention=np.nanmedian(v) if len(v) else np.nan,
                                 p_sign_vs_1=p,
                                 match_scale_um=np.nanmax(
                                     [match_scale(markers.get(dep, {}).get(d),
                                                  TOOLS[t][slot])
                                      for d in sorted(markers.get(dep, {}))] or [np.nan])))
        for i in range(len(names)):
            for arr, off in ((lo_r, off_lo), (hi_r, off_hi)):
                if not arr.size or not np.isfinite(arr[:, i]).any():
                    ax.text(i + off, 0.0, "n.d.", ha="center", va="bottom",
                            fontsize=5.4, color=INK2, rotation=90, zorder=8)
        if j == 0:
            ax.set_ylabel("Marker retention ratio\n(log$_2$ scale)", fontweight="bold",
                          labelpad=4)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=45, ha="right", rotation_mode="anchor",
                           fontsize=6.2)
        axs[0, j].set_xticks(x); axs[0, j].set_xticklabels([])
        for a_ in (axs[0, j], axs[1, j]):
            a_.set_xlim(-0.55, len(names) - 0.45)
            a_.yaxis.grid(True, color=GRID, lw=0.5, zorder=0); a_.set_axisbelow(True)
    for j in range(len(DEPTHS)):
        axs[0, j].set_ylim(0, max(top_a * 1.15, 1.0))
        pad = 0.05 * (ret_hi - ret_lo)
        axs[1, j].set_ylim(ret_lo - pad, ret_hi + pad)
        tk = [t for t in range(-4, 5) if ret_lo - pad <= t <= ret_hi + pad]
        axs[1, j].set_yticks(tk)
        axs[1, j].set_yticklabels([f"{2.0 ** t:g}" for t in tk])
    if key:
        fig.legend(handles=[Patch(facecolor=C_ALONE, label=KEY_NATIVE),
                            Patch(facecolor=C_DECONV, label=KEY_AD)],
                   loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0),
                   handlelength=1.3, columnspacing=1.6)
    fig.tight_layout(w_pad=1.4, h_pad=1.0, rect=(0, 0.06 if key else 0, 1, 1))
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{out}.{ext}", bbox_inches="tight")
        print("wrote", f"{out}.{ext}")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(f"{out}_stats.csv", index=False)

def depth_trend(scores, out):
    names = list(TOOLS)
    x = np.arange(len(DEPTHS))
    fig, axs = plt.subplots(len(P.CHEMS), len(names), figsize=(7.2, 4.9), sharey=True)
    top = max((discard_pct(scores[d][ds], TOOLS[t][sl])
               for d in DEPTHS for ds in scores.get(d, {})
               for t in names for sl in (0, 1)), default=1.0)
    nper = {}
    for r, ch in enumerate(P.CHEMS):
        for i, t in enumerate(names):
            ax = axs[r, i]
            ends = []
            for slot, col, lab in ((0, C_ALONE, KEY_NATIVE), (1, C_DECONV, KEY_AD)):
                per = [[discard_pct(scores[d][ds], TOOLS[t][slot])
                        for ds in sorted(scores.get(d, {})) if P.chemistry(ds) == ch]
                       for d in DEPTHS]
                nper[ch] = max(len(v) for v in per)
                m = np.array([np.nanmedian(v) if v else np.nan for v in per], float)
                for j, v in enumerate(per):
                    if v:
                        jit = np.linspace(-.09, .09, len(v))
                        ax.scatter(x[j] + jit, v, s=5, color="#ffffff", edgecolors=col,
                                   linewidths=.7, zorder=4)
                ax.plot(x, m, color=col, lw=2, marker="o", ms=4.5, mec="white", mew=.8,
                        zorder=5, label=lab)
                if np.isfinite(m[-1]):
                    ends.append((m[-1], col))
            if len(ends) == 2 and abs(ends[0][0] - ends[1][0]) < 0.06 * max(top, 1.0):
                ends = [(ends[0][0], ends[0][1], 7), (ends[1][0], ends[1][1], -7)]
            else:
                ends = [(v, c, 0) for v, c in ends]
            for v, col, dy in ends:
                ax.annotate(f"{v:.0f}%", (x[-1], v), xytext=(6, dy),
                            textcoords="offset points", va="center", fontsize=7,
                            color=col, fontweight="bold")
            if r == 0:
                ax.set_title(t, pad=5, color=INK)
            ax.set_xticks(x)
            ax.set_xticklabels([str(d) for d in DEPTHS])
            ax.set_xlim(-0.35, len(DEPTHS) - 0.35)
            if r == len(P.CHEMS) - 1:
                ax.set_xlabel("UMIs per bin", fontsize=6.8, color=INK2)
            ax.yaxis.grid(True, color=GRID, lw=.5, zorder=0); ax.set_axisbelow(True)
        axs[r, 0].set_ylabel("Undamaged tissue\nremoved (%)")
        axs[r, -1].text(1.05, 0.5, P.CHEM_LABEL[ch], transform=axs[r, -1].transAxes,
                        rotation=270, ha="left", va="center", fontsize=7.5,
                        fontweight="bold", color=INK)
    axs[0, 0].set_ylim(0, top * 1.12)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=False, handlelength=1.6,
               columnspacing=2.4, fontsize=7, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(w_pad=1.6, h_pad=1.4, rect=(0, 0, 1, 0.95))
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{out}.{ext}", bbox_inches="tight")
        print("wrote", f"{out}.{ext}")
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="/work/benchmark_grid")
    ap.add_argument("--outdir", default="/work/benchmark_grid/figures")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    style()
    scores, markers = load_grid(a.grid)
    for ch in P.CHEMS:
        keep = lambda tbl: {d: {k: v for k, v in tbl[d].items() if P.chemistry(k) == ch}
                            for d in tbl}
        figure(keep(scores), keep(markers),
               os.path.join(a.outdir, f"fig2_discarded_tissue_{ch}"), key=(ch == P.CHEMS[-1]))
    depth_trend(scores, os.path.join(a.outdir, "fig3_depth_trend"))

if __name__ == "__main__":
    main()
