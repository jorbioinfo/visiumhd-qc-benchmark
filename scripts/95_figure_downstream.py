#!/usr/bin/env python3
import argparse, importlib.util, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import transforms
from matplotlib.colors import ListedColormap

import panel as P
import figstyle as FS

REF = ("01_no_qc", "No filter", "#6b7280")
ARMS = [("02_mad_8um", "native\n8 µm", "#2a78d6"),
        ("10_mad_64um", "+ A&P\n64 µm", "#eb6834"),
        (P.SPOTSWEEPER[8], "native\n8 µm", "#2a78d6"),
        (P.SPOTSWEEPER[64], "+ A&P\n64 µm", "#eb6834"),
        (P.SPATIALARTIFACTS[8], "native\n8 µm", "#2a78d6"),
        (P.SPATIALARTIFACTS[16], "+ A&P\n16 µm", "#eb6834")]
GROUPS = [("MAD", 0, 1), ("SpotSweeper", 2, 3), ("SpatialArtifacts", 4, 5)]
W = 7.2

MAP_SECTION = "Visium_HD_Human_Colon_Cancer"
MAP_ARMS = [("01_no_qc", "No filter"), ("02_mad_8um", "MAD, native 8 µm bins"),
            ("10_mad_64um", "MAD + A&P, 64 µm bins")]
MAP_DEPTHS = (100, 12)
INK = "#1f2a24"

def maps_module():
    spec = importlib.util.spec_from_file_location(
        "maps122", os.path.join(os.path.dirname(os.path.abspath(__file__)), "122_demo_downstream_maps.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def maps(fig, cell, labels, summary, grid):
    M = maps_module()
    df = pd.read_csv(labels, dtype={"label": str})
    sm = pd.read_csv(summary).query("section == @MAP_SECTION").set_index("arm")
    rem, _frac = M.removal(MAP_SECTION, grid)
    allsq = pd.Index(sorted(set(df[df.arm == REF[0]].square)))
    r, c = M.rc_of(allsq)
    R, C = r.max() + 1, c.max() + 1
    pos = {sq: (i, j) for sq, i, j in zip(allsq, r, c)}
    pal = M.palette(df["label"].unique())
    purples = plt.get_cmap("Purples")
    nrow = len(MAP_DEPTHS) + 1
    gs = cell.subgridspec(nrow, len(MAP_ARMS), wspace=0.04, hspace=0.26)
    axs = np.array([[fig.add_subplot(gs[i, j]) for j in range(len(MAP_ARMS))]
                    for i in range(nrow)])

    def img(values, colour):
        im = np.full((R, C, 4), 1.0)
        for sq in allsq:
            im[pos[sq]] = (0.80, 0.80, 0.80, 1)
        for sq, v in values.items():
            if sq in pos:
                im[pos[sq]] = colour(v)
        return im

    rows = []
    for j, (arm, name) in enumerate(MAP_ARMS):
        for i, d in enumerate(MAP_DEPTHS):
            s = df[(df.arm == arm) & (df.depth == d)].set_index("square")["label"]
            ax = axs[i, j]
            ax.imshow(img(s, lambda v: pal[v]), interpolation="nearest")
            removed = 100 * sm.loc[arm, f"removed_{d}"]
            ax.set_title(f"{s.nunique()} clusters, {removed:.1f}% removed", fontsize=6.2,
                         color=INK, pad=2.5)
            rows.append(dict(section=MAP_SECTION, arm=arm, depth=d, clusters=s.nunique(),
                             removed_pct=round(removed, 1), ari=sm.loc[arm, "ari"]))
        axs[0, j].text(0.5, 1.2, name, transform=axs[0, j].transAxes, ha="center", va="bottom",
                       fontsize=7, fontweight="bold", color=INK)
        axs[1, j].set_xlabel(f"ARI {sm.loc[arm, 'ari']:.3f}", fontsize=6.8, color=INK, labelpad=3)
        axs[2, j].imshow(img(rem[arm, MAP_DEPTHS[-1]], lambda v: purples(0.12 + 0.88 * v)),
                         interpolation="nearest")
    for ax in axs.ravel():
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    for i, d in enumerate(MAP_DEPTHS):
        axs[i, 0].set_ylabel(f"{d} UMIs\nper bin", fontsize=7, fontweight="bold", rotation=0,
                             ha="right", va="center", labelpad=6, color=INK, linespacing=1.2)
    axs[2, 0].set_ylabel(f"Removed at\n{MAP_DEPTHS[-1]} UMIs\nper bin", fontsize=7, fontweight="bold",
                         rotation=0, ha="right", va="center", labelpad=6, color=INK, linespacing=1.2)
    cax = axs[2, 1].inset_axes([0.1, -0.13, 0.8, 0.05])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 100), cmap=ListedColormap(
        purples(np.linspace(0.12, 1, 256)))), cax=cax, orientation="horizontal")
    cb.set_ticks([0, 50, 100]); cb.set_ticklabels(["0%", "50%", "100%"])
    cb.set_label("share of the 32 µm square removed", fontsize=5.8, color=INK, labelpad=1.5)
    cb.ax.tick_params(labelsize=5.6, colors=INK, length=1.5, pad=1); cb.outline.set_visible(False)
    return axs, pd.DataFrame(rows)

def group_labels(ax, fig, off):
    base = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    line = base + transforms.ScaledTranslation(0, -25 / 72, fig.dpi_scale_trans)
    text = base + transforms.ScaledTranslation(0, -27 / 72, fig.dpi_scale_trans)
    for name, i0, i1 in GROUPS:
        ax.plot([off + i0 - 0.34, off + i1 + 0.34], [0, 0], transform=line,
                color="#555555", lw=0.6, clip_on=False)
        ax.text(off + (i0 + i1) / 2, 0, name, transform=text, ha="center", va="top",
                fontsize=6.8)

def spread(i, n):
    return np.full(n, i) + (np.linspace(-0.15, 0.15, n) if n > 1 else [0.0])

def bars(ax, i, v, clr):
    ax.bar(i, np.median(v), width=0.62, color=clr, alpha=0.22,
           edgecolor=clr, linewidth=1.0, zorder=2)
    ax.scatter(spread(i, len(v)), v, s=9, color=clr, zorder=4,
               edgecolor="white", linewidth=0.4)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="benchmark_grid/figures/downstream_binned.csv")
    ap.add_argument("--labels", default="benchmark_grid/figures/demo_downstream/"
                                        f"{MAP_SECTION}_labels.csv.gz",
                    help="labels 122_demo_downstream_maps.py saved for the section of panel a")
    ap.add_argument("--summary", default="benchmark_grid/figures/demo_downstream/summary.csv")
    ap.add_argument("--grid", default="benchmark_grid",
                    help="the benchmark grid, for the share of each square removed in panel a")
    ap.add_argument("--out", default="benchmark_grid/figures/fig4_clustering_agreement")
    a = ap.parse_args()
    d = pd.read_csv(a.csv).drop_duplicates(
        subset=["section", "arm", "comparison"], keep="last")
    d = d[d.section.map(lambda x: P.short(x) in P.SECTIONS)]
    d["chem"] = d.section.map(P.chemistry)
    cols = [REF] + ARMS
    keep = [c[0] for c in cols]

    plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.3,
                         "ytick.labelsize": 6.5, "legend.fontsize": 6.8,
                         "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(W, 9.3))
    outer = fig.add_gridspec(2, 1, height_ratios=[4.9, 4.4], hspace=0.14)
    map_axs, map_stats = maps(fig, outer[0], a.labels, a.summary, a.grid)
    inner = outer[1].subgridspec(len(P.CHEMS), 2, width_ratios=[len(cols), len(ARMS)],
                                 wspace=0.32, hspace=0.12)
    ax = np.array([[fig.add_subplot(inner[r, k]) for k in range(2)] for r in range(len(P.CHEMS))])
    for k in range(2):
        for r in range(1, len(P.CHEMS)):
            ax[r, k].sharey(ax[0, k])
    out, ns = [], {}
    for r, ch in enumerate(P.CHEMS):
        dc = d[d.chem == ch]
        ad = dc[dc.comparison == "across_depth"].pivot_table(
            index="section", columns="arm", values="ari")
        rs = dc[dc.comparison == "reseed"].pivot_table(
            index="section", columns="arm", values="ari")
        ad = ad[[k for k in keep if k in ad.columns]].dropna()
        ns[ch] = len(ad)
        A, B = ax[r, 0], ax[r, 1]
        for i, (arm, lab, clr) in enumerate(cols):
            if arm in ad.columns:
                bars(A, i, ad[arm].to_numpy(), clr)
        for i, (arm, lab, clr) in enumerate(ARMS):
            if arm in ad.columns:
                bars(B, i, (ad[arm] - ad[REF[0]]).to_numpy(), clr)
        A.axhline(0, lw=0.8, color="#999999", zorder=1)
        B.axhline(0, ls="--", lw=1.0, color="#555555", zorder=3)
        A.set_ylabel("ARI between the two depths")
        B.set_ylabel("ARI difference from no filter")
        B.text(1.04, 0.5, P.CHEM_LABEL[ch], transform=B.transAxes, rotation=270,
               ha="left", va="center", fontsize=7.5, fontweight="bold")
        A.set_xticks(range(len(cols))); B.set_xticks(range(len(ARMS)))
        last = r == len(P.CHEMS) - 1
        A.set_xticklabels([c[1] for c in cols] if last else [])
        B.set_xticklabels([c[1] for c in ARMS] if last else [])
        if last:
            group_labels(A, fig, 1)
            group_labels(B, fig, 0)
        t = ad.copy(); t.insert(0, "chem", ch)
        for arm in [c[0] for c in cols]:
            if arm in rs.columns:
                t[f"reseed_{arm}"] = rs.loc[ad.index, arm]
        out.append(t)
    for x in ax.ravel():
        x.grid(axis="y", lw=0.5, color="#dddddd", zorder=0)
        x.set_axisbelow(True)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.tight_layout(h_pad=1.2, w_pad=2.5)
    FS.letter_axes(fig, map_axs[0, 0], "a")
    FS.letter_axes(fig, ax[0, 0], "b")
    FS.letter_axes(fig, ax[0, 1], "c")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{a.out}.{ext}", dpi=600, bbox_inches="tight")

    map_stats.to_csv(a.out + "_maps.csv", index=False)
    print("panel a:")
    print(map_stats.to_string(index=False))
    res = pd.concat(out)
    res.to_csv(a.out + "_stats.csv")
    print(res.round(3).to_string())
    print("\ndifference from no filter:")
    arms = [c[0] for c in ARMS if c[0] in res.columns]
    diff = res[arms].sub(res[REF[0]], axis=0)
    diff.insert(0, "chem", res.chem)
    print(diff.round(3).to_string())

if __name__ == "__main__":
    main()
