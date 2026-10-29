import argparse, glob, os, re

import panel as P
import numpy as np, pandas as pd
from scipy.io import mmread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

INK, INK2, MUTED, BASE = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
HUE = "#2a78d6"
CAT = ["#2a78d6", "#eb6834", "#1baf7a"]

def style():
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
        "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8.5,
        "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "legend.fontsize": 6.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": BASE, "axes.linewidth": 0.6,
        "figure.facecolor": "white", "axes.facecolor": "white"})

HUE2 = "#eb6834"

def mito_baseline(src):
    rows = []
    for sec in sorted(P.SECTIONS):
        name = P.PREFIX + sec
        d = os.path.join(src, name, "binned_outputs", "square_008um",
                         "filtered_feature_bc_matrix")
        feat = pd.read_csv(os.path.join(d, "features.tsv.gz"), sep="\t", header=None)
        sym = feat[1].astype(str)
        mt = sym.str.match(r"^(MT-|mt-)", na=False).to_numpy()
        if not mt.any():
            print("no MT genes in", name, flush=True); continue
        m = mmread(os.path.join(d, "matrix.mtx.gz")).tocsr()
        tot = np.asarray(m.sum(0)).ravel().astype(float)
        mito = np.asarray(m[mt].sum(0)).ravel().astype(float)
        keep = tot > 0
        frac = mito[keep] / tot[keep] * 100
        rows.append(dict(dataset=sec.replace("_", " "), chem=P.chemistry(sec),
                         median=float(np.median(frac)),
                         p90=float(np.percentile(frac, 90)), n=int(keep.sum())))
        print(f"{name:<40} median {rows[-1]['median']:.1f}%  p90 {rows[-1]['p90']:.1f}%",
              flush=True)
    return order_s4(pd.DataFrame(rows))

def order_s4(df):
    rank = {c: i for i, c in enumerate(P.CHEMS)}
    return df.assign(_r=df.chem.map(rank)).sort_values(["_r", "median"]).drop(columns="_r")

def fig_s4(df, outdir):
    fig, ax = plt.subplots(figsize=(5.2, 0.17 * len(df) + 1.4))
    y = np.arange(len(df)) + (df.chem.to_numpy() == P.POLYA)
    col = {P.PROBE: HUE, P.POLYA: HUE2}
    ax.barh(y, df["median"], height=0.62, color=[col[c] for c in df.chem], zorder=3)

    ax.plot(df["p90"], y, "|", ms=7, mew=1.3, color=INK2, zorder=4)
    for cut, lab in ((10, "10%"), (20, "20%")):
        ax.axvline(cut, color=BASE, lw=0.8, ls=(0, (3, 2)), zorder=2)
        ax.text(cut, y.max() + 0.65, lab, ha="center", va="bottom",
                fontsize=6.2, color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels([d.replace("3prime ", "3' ").replace("11mm ", "11 mm ")
                        for d in df["dataset"]], fontsize=6.8)
    ax.set_xlabel("Mitochondrial fraction per 8 µm bin (%)", color=INK)
    ax.tick_params(colors=INK2)
    ax.legend([Patch(color=col[P.PROBE]), Patch(color=col[P.POLYA]),
               Line2D([], [], marker="|", ls="none", ms=7, mew=1.3, color=INK2)],
              [f"median, {P.CHEM_LABEL[P.PROBE]}", f"median, {P.CHEM_LABEL[P.POLYA]}",
               "90th percentile"], frameon=False, loc="upper center", ncol=3,
              bbox_to_anchor=(0.5, -0.12), handlelength=1.2)
    fig.tight_layout()
    for e in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(outdir, f"figS4_mito_baseline.{e}"), bbox_inches="tight")
    plt.close(fig)

ARM_LABEL = {
    "02_mad_8um": "MAD, native 8 µm", "08_mad_16um": "MAD + A&D, 16 µm",
    "09_mad_32um": "MAD + A&D, 32 µm", "10_mad_64um": "MAD + A&D, 64 µm",
    "41_mad_128um": "MAD + A&D, 128 µm",
    "35_mad_64um_excludepartial": "MAD + A&D, 64 µm, partly covered bins left out",
    "24_spotsweeper_8um_k48": "SpotSweeper, native 8 µm", "25_spotsweeper_16um_k48": "SpotSweeper + A&D, 16 µm",
    "26_spotsweeper_32um_k48": "SpotSweeper + A&D, 32 µm", "27_spotsweeper_64um_k48": "SpotSweeper + A&D, 64 µm",
    "28_spotsweeper_128um_k48": "SpotSweeper + A&D, 128 µm",
    "44_spotsweeper_8um_k192": "SpotSweeper, native 8 µm, 192 neighbours",
    "45_spotsweeper_8um_k768": "SpotSweeper, native 8 µm, 768 neighbours",
    "46_spotsweeper_8um_k3072": "SpotSweeper, native 8 µm, 3072 neighbours",
    "03_spotsweeper_8um": "SpotSweeper, native 8 µm, 36 neighbours",
    "30_miqc_8um": "miQC, native 8 µm", "31_miqc_64um": "miQC + A&D, 64 µm",
    "19_spatialartifacts_8um": "SpatialArtifacts, native 8 µm",
    "20_spatialartifacts_16um": "SpatialArtifacts + A&D, 16 µm",
    "32_fixedumi_8um": "Fixed count, native 8 µm", "33_fixedumi_64um": "Fixed count + A&D, 64 µm",
}

def fig_s6(grid, outdir):
    rows = []
    for _sec, _dep, cell in P.cells(grid, need=("scored", "marker.csv")):
        p = os.path.join(cell, "scored", "marker.csv")
        f = pd.read_csv(p)
        f["run"] = os.path.basename(cell)
        rows.append(f)
    d = pd.concat(rows)
    print("S6 arms left out:", sorted(set(d.arm) - set(ARM_LABEL)))
    d = d[d.arm.isin(ARM_LABEL)]
    lvl = {(256.0, 20.0): "256 µm, 20 strata", (512.0, 10.0): "512 µm, 10 strata",
           (1024.0, 5.0): "1024 µm, 5 strata"}
    d["level"] = ["unmatched" if not np.isfinite(a) else lvl.get((a, b), "other")
                  for a, b in zip(d.match_scale_um, d.match_strata)]
    assert not (d.level == "other").any(), d[d.level == "other"].head()
    order = list(lvl.values()) + ["unmatched"]
    arms = [a for a in ARM_LABEL if a in set(d.arm)][::-1]
    counts = pd.DataFrame({k: [((d.arm == a) & (d.level == k)).sum() for a in arms]
                           for k in order}, index=arms)
    fig, ax = plt.subplots(figsize=(6.0, max(2.4, 0.26 * len(arms) + 1.2)))
    left = np.zeros(len(arms))
    for k, c in zip(order, CAT + ["#c3c2b7"]):
        ax.barh(np.arange(len(arms)), counts[k], left=left, height=0.62,
                color=c, label=k, zorder=3)
        left = left + counts[k].to_numpy()
    ax.set_yticks(np.arange(len(arms)))
    ax.set_yticklabels([ARM_LABEL[a] for a in arms], fontsize=6.2)
    ax.set_xlabel("Runs", color=INK)
    ax.tick_params(colors=INK2)
    ax.legend(frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.08), title="Matching level", handlelength=1.2,
              columnspacing=1.2)
    fig.tight_layout()
    for e in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(outdir, f"figS6_matching_ladder.{e}"), bbox_inches="tight")
    plt.close(fig)
    print("S6 totals:", {k: int(counts[k].sum()) for k in order})

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="/work/public")
    ap.add_argument("--grid", default="/work/benchmark_grid")
    ap.add_argument("--outdir", default="/work/benchmark_grid/figures")
    ap.add_argument("--from-csv", action="store_true",
                    help="redraw S4 from the saved table instead of re-reading matrices")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    style()
    csv = os.path.join(a.outdir, "figS4_mito_baseline.csv")
    fig_s6(a.grid, a.outdir)
    print("wrote figS6_matching_ladder", flush=True)
    if a.from_csv:
        df = pd.read_csv(csv)
        df = order_s4(df[df.dataset.str.replace(" ", "_").isin(P.SECTIONS)])
    else:
        df = mito_baseline(a.src)
        df.to_csv(csv, index=False)
    fig_s4(df, a.outdir)
    print("wrote figS4_mito_baseline")

if __name__ == "__main__":
    main()
