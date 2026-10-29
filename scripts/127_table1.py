import csv, os

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(WK, "benchmark_grid/figures/table1_panel.csv")

PREP = {
    "11mm_Human_Breast_Cancer": "not stated", "11mm_Human_Colon_Cancer_Fresh_Frozen": "fresh frozen",
    "11mm_Human_Pancreas": "not stated", "11mm_Mouse_Embryo": "not stated",
    "6p5mm_Human_Heart": "FFPE", "Human_Breast_Cancer_FFPE": "FFPE",
    "Human_Breast_Cancer_Fresh_Frozen": "fresh frozen", "Human_Colon_Cancer": "FFPE",
    "Human_Colon_Cancer_P1": "FFPE", "Human_Colon_Cancer_P5": "FFPE",
    "Human_Colon_Normal_P3": "FFPE", "Human_Colon_Normal_P5": "FFPE",
    "Human_Kidney_FFPE": "FFPE", "Human_Lung_Cancer": "FFPE",
    "Human_Lung_Cancer_Fixed_Frozen": "fixed frozen", "Human_Lymph_Node_FFPE": "FFPE",
    "Human_Ovarian_Cancer_FF": "fresh frozen", "Human_Pancreas": "FFPE",
    "Human_Prostate_Cancer_FFPE": "FFPE", "Human_Tonsil_Fresh_Frozen": "fresh frozen",
    "Human_Tonsil_Fresh_Frozen_IF": "fresh frozen", "Mouse_Brain": "FFPE",
    "Mouse_Brain_Fixed_Frozen": "fixed frozen", "Mouse_Brain_Fresh_Frozen": "fresh frozen",
    "Mouse_Embryo": "FFPE", "Mouse_Kidney": "FFPE", "Mouse_Lung_Fresh_Frozen": "fresh frozen",
    "Mouse_Small_Intestine": "FFPE",
    "SixCancer_Bladder_Cancer": "FFPE", "SixCancer_Colorectal_Cancer": "FFPE",
    "SixCancer_Lung_Cancer": "FFPE", "SixCancer_Lymphoma": "FFPE", "SixCancer_Ovarian_Cancer": "FFPE",
}
SPACERANGER = {
    "11mm_Human_Pancreas": "4.1.0", "Mouse_Brain_Fresh_Frozen": "3.1.1", "11mm_Mouse_Embryo": "4.1.0",
    "Human_Breast_Cancer_Fresh_Frozen": "3.1.1", "Mouse_Embryo": "4.0.1", "Mouse_Lung_Fresh_Frozen": "3.1.1",
    "11mm_Human_Colon_Cancer_Fresh_Frozen": "4.1.0", "Human_Pancreas": "4.0.1",
    "Human_Tonsil_Fresh_Frozen": "3.1.1", "Human_Colon_Cancer": "4.0.1", "Mouse_Kidney": "4.0.1",
    "11mm_Human_Breast_Cancer": "4.1.0", "Mouse_Small_Intestine": "4.0.1", "Human_Lung_Cancer": "3.0.0",
    "Mouse_Brain_Fixed_Frozen": "3.1.1", "Human_Tonsil_Fresh_Frozen_IF": "3.1.1",
    "Human_Lung_Cancer_Fixed_Frozen": "3.1.1", "Mouse_Brain": "4.0.1", "Human_Kidney_FFPE": "4.0.1",
    "Human_Lymph_Node_FFPE": "4.0.1", "6p5mm_Human_Heart": "4.1.0", "Human_Prostate_Cancer_FFPE": "4.0.1",
    "3prime_Mouse_Kidney": "4.0.1", "3prime_Mouse_Ovary": "4.1.0", "3prime_Mouse_Embryo": "4.0.1",
    "3prime_Human_Tonsil": "4.1.0", "3prime_Mouse_Brain": "4.0.1", "3prime_Human_Lymph_Node": "4.1.0",
    "3prime_Human_Ovarian_Cancer": "4.0.1", "Human_Breast_Cancer_FFPE": "3.1.2",
    "Human_Ovarian_Cancer_FF": "4.0.1",
    **{s: "3.0.0" for s in ("Human_Colon_Cancer_P1", "Human_Colon_Cancer_P5", "Human_Colon_Normal_P3",
                            "Human_Colon_Normal_P5", "SixCancer_Bladder_Cancer",
                            "SixCancer_Colorectal_Cancer", "SixCancer_Lung_Cancer",
                            "SixCancer_Lymphoma", "SixCancer_Ovarian_Cancer")},
}
SOURCE = {"Human_Colon_Cancer_P1": "[3]", "Human_Colon_Cancer_P5": "[3]",
          "Human_Colon_Normal_P3": "[3]", "Human_Colon_Normal_P5": "[3]",
          **{f"SixCancer_{t}": "[16]" for t in ("Bladder_Cancer", "Colorectal_Cancer", "Lung_Cancer",
                                               "Lymphoma", "Ovarian_Cancer")}}

def label(sec):
    s = sec
    for a, b in (("11mm_", ""), ("6p5mm_", ""), ("3prime_", ""), ("SixCancer_", "Six-cancer "),
                 ("_Fresh_Frozen_IF", ", immunofluorescence"), ("_Fresh_Frozen", ""),
                 ("_Fixed_Frozen", ""), ("_FFPE", ""), ("_FF", "")):
        s = s.replace(a, b)
    return s.replace("_", " ")

def prep(r):
    return "fresh frozen" if r["chemistry"] == "3prime" else PREP.get(r["section"], "not stated")

if __name__ == "__main__":
    rows = list(csv.DictReader(open(SRC)))
    for r in rows:
        r["prep"] = prep(r)
    rows.sort(key=lambda r: (r["chemistry"] != "probe", -float(r["native_mean_umi"] or 0)))
    print("| Section | Species | Chemistry | Preparation | Slide | Space Ranger | In-tissue 8 µm bins "
          "| Native UMIs/bin | Mito. baseline |")
    print("|---|---|---|---|---|---|---:|---:|---:|")
    for r in rows:
        chem = "probe" if r["chemistry"] == "probe" else "3′ polyA"
        slide = "11 mm" if float(r["slide_mm"]) == 11.0 else "6.5 mm"
        src = SOURCE.get(r["section"], "")
        name = label(r["section"]) + (f" {src}" if src else "")
        bins = f"{int(r['in_tissue_bins']):,}" if r["in_tissue_bins"] else "-"
        umi = f"{float(r['native_mean_umi']):,.0f}" if r["native_mean_umi"] else "-"
        mt = f"{float(r['mito_baseline_pct']):.1f}%" if r["mito_baseline_pct"] else "-"
        print(f"| {name} | {r['species']} | {chem} | {r['prep']} | {slide} | {r['spaceranger'] or SPACERANGER.get(r['section'], '-')} "
              f"| {bins} | {umi} | {mt} |")
