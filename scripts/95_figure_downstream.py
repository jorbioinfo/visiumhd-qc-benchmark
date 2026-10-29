import argparse, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import panel as P
import figstyle as FS

REF = ("01_no_qc", "No filter", "#6b7280")
ARMS = [("02_mad_8um", "MAD\nnative 8 µm", "#2a78d6"),
        ("10_mad_64um", "MAD + A&D\n64 µm bins", "#eb6834"),
        (P.SPOTSWEEPER[8], "SpotSweeper\nnative 8 µm", "#1baf7a")]
W = 7.2

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
    ap.add_argument("--out", default="benchmark_grid/figures/fig8_downstream")
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
    fig, ax = plt.subplots(len(P.CHEMS), 2, figsize=(W, 5.4), sharey="col")
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
                if arm in rs.columns:
                    A.scatter([i], [rs.loc[ad.index, arm].median()], marker="D", s=16,
                              facecolor="none", edgecolor="#333333", linewidth=0.9,
                              zorder=5)
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
        t = ad.copy(); t.insert(0, "chem", ch); out.append(t)
    ax[-1, 0].scatter([], [], marker="D", s=16, facecolor="none", edgecolor="#333333",
                      linewidth=0.9, label="same survivors, second seed (median)")
    ax[-1, 0].legend(frameon=False, loc="upper center",
                     bbox_to_anchor=(0.5, -0.2), handletextpad=0.4)
    for x in ax.ravel():
        x.grid(axis="y", lw=0.5, color="#dddddd", zorder=0)
        x.set_axisbelow(True)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.tight_layout(w_pad=2.5)
    FS.letter_axes(fig, ax[0, 0], "a")
    FS.letter_axes(fig, ax[0, 1], "b")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{a.out}.{ext}", dpi=600, bbox_inches="tight")

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
