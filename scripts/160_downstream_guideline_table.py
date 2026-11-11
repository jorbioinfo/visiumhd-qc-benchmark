#!/usr/bin/env python3
import math, os, sys
from decimal import Decimal, ROUND_HALF_UP
import pandas as pd

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(WK, "scripts"))
import panel as P

FIG = os.path.join(WK, "benchmark_grid", "figures")
ARM = {"none": "none", "02_mad_8um": "MAD, 8 µm", "09_mad_32um": "MAD + A&P, 32 µm",
       "10_mad_64um": "MAD + A&P, 64 µm",
       "44_spotsweeper_8um_k192": "SpotSweeper, 8 µm, 192 nb.",
       "45_spotsweeper_8um_k768": "SpotSweeper, 8 µm, 768 nb.",
       "25_spotsweeper_16um_k48": "SpotSweeper + A&P, 16 µm",
       "26_spotsweeper_32um_k48": "SpotSweeper + A&P, 32 µm",
       "27_spotsweeper_64um_k48": "SpotSweeper + A&P, 64 µm",
       "28_spotsweeper_128um_k48": "SpotSweeper + A&P, 128 µm",
       "19_spatialartifacts_8um": "SpatialArtifacts, 8 µm",
       "20_spatialartifacts_16um": "SpatialArtifacts + A&P, 16 µm"}
ROWS = [("no_filter", "No filter"), ("mad8_both", "MAD, native bins, at both depths"),
        ("mad_64um_both", "MAD, 64 µm, at both depths"), ("procedure_section_wide", "Section-wide test"),
        ("procedure_0_45um", "Damage < 45 µm"), ("procedure_45_130um", "Damage 45–130 µm"),
        ("procedure_130_360um", "Damage 130–360 µm"), ("procedure_360_upum", "Damage > 360 µm")]

def sign_p(up, dn):
    n = up + dn
    return float("nan") if n == 0 else min(1.0, 2 * sum(math.comb(n, k) for k in range(min(up, dn) + 1)) / 2 ** n)

def fp(p):
    if p < 0.001:
        return "< 0.001"
    if p >= 0.995:
        return "1.0"
    d = Decimal(repr(p))
    return str(d.quantize(Decimal(1).scaleb(d.adjusted() - 1), rounding=ROUND_HALF_UP))

def label(arms):
    return "; ".join(ARM[a] for a in arms.split("+"))

G = pd.read_csv(os.path.join(FIG, "downstream_guideline.csv"))
B = pd.read_csv(os.path.join(FIG, "downstream_binned.csv"))
g = G[G.comparison == "across_depth"]
ari = g.pivot(index="section", columns="arm", values="ari")
b = B[B.comparison == "across_depth"]
ari["mad8_both"] = b.pivot(index="section", columns="arm", values="ari")["02_mad_8um"]
arms = {a: (label(x.arms_a.iloc[0]), label(x.arms_b.iloc[0])) for a, x in g.groupby("arm")}
for a, x in g.groupby("arm"):
    assert x.arms_a.nunique() == 1 and x.arms_b.nunique() == 1, a
arms["mad8_both"] = (ARM["02_mad_8um"], ARM["02_mad_8um"])
ari["chem"] = [P.chemistry(s) for s in ari.index]

def versus(x, a, ref, lower):
    y = x[[a, ref]].dropna()
    k = int((y[a] < y[ref]).sum()) if lower else int((y[a] > y[ref]).sum())
    m = int((y[a] > y[ref]).sum()) if lower else int((y[a] < y[ref]).sum())
    return f"{k} of {len(y)} ({fp(sign_p(k, m))})"

out = ["| Chemistry | Setting | 100 UMIs/bin | 12 UMIs/bin | n | ARI | Change from no filter | "
       "Lower than no filter | Higher than MAD, native bins | Lower than MAD, 64 µm |",
       "|---|---|---|---|---:|---:|---:|---:|---:|---:|"]
for chem, cname in (("probe", "Probe-based"), ("3prime", "3′ polyA")):
    x = ari[ari.chem == chem]
    first = True
    for a, name in ROWS:
        v = x[a].dropna()
        d = (x[a] - x["no_filter"]).dropna()
        cells = [cname if first else "", name, arms[a][0] if a != "no_filter" else "–",
                 arms[a][1] if a != "no_filter" else "–", str(len(v)), f"{v.median():.3f}"]
        if a == "no_filter":
            cells += ["–", "–", versus(x, a, "mad8_both", False), versus(x, a, "mad_64um_both", True)]
        else:
            cells += [f"{d.median():+.3f}".replace("-", "−"), versus(x, a, "no_filter", True),
                      "–" if a == "mad8_both" else versus(x, a, "mad8_both", False),
                      "–" if a == "mad_64um_both" else versus(x, a, "mad_64um_both", True)]
        out.append("| " + " | ".join(cells) + " |")
        first = False
md = "\n".join(out)
open(os.path.join(FIG, "downstream_guideline_table.md"), "w").write(md + "\n")
print(md)
