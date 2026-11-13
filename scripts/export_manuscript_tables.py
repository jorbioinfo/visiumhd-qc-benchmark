#!/usr/bin/env python3
"""Export the manuscript tables and optionally package its final figures and CSVs."""
import argparse
import csv
from pathlib import Path
import re
from shutil import copyfile
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "benchmark_grid" / "figures"
CAPTION = re.compile(r"^\*\*Table (S?\d+)\.")
FIGURES = (
    "fig1_study_design", "fig2_binsize_recall", "fig3_depth_trend",
    "fig4_clustering_agreement", "fig5_binsize_f1", "fig6_qc_guideline",
    "fig7_qc_arms", "fig8_synthetic_lqrs", "fig9_clustering_design",
    "figS1_recall_mito", "figS2_recall_signal_loss", "figS3_marker_retention",
    "figS4_false_positive_rate", *(f"figS5{panel}_clustering_maps" for panel in "abcd"),
    "figS6a_mad", "figS6b_miqc", "figS6c_spotsweeper", "figS6d_spatialartifacts",
    "figS7_binsize_d50", "figS8_binsize_d25", "figS9_molecule_curve",
    "figS10_widened_neighbourhood", "figS11_mito_baseline", "figS12_marker_check",
    "figS13_matching_ladder",
)
DIAGRAMS = {"fig1_study_design", "fig7_qc_arms", "fig8_synthetic_lqrs",
            "fig9_clustering_design", "figS12_marker_check"}


def manuscript_tables(drafts):
    tables = {}
    for draft in drafts:
        blocks = draft.read_text(encoding="utf-8").split("\n\n")
        for index, (caption, body) in enumerate(zip(blocks, blocks[1:])):
            match = CAPTION.match(caption)
            if match and body.startswith("|"):
                number = match[1]
                assert number not in tables, number
                note = ""
                if number == "1" and blocks[index + 2].startswith("Values are "):
                    note = f"\n{blocks[index + 2]}\n"
                tables[number] = (caption, body, note)

    # Table S1 is generated separately by 82_section_level_tests.py.
    assert set(tables) == {"1", *(f"S{i}" for i in range(2, 12))}, tables.keys()
    return tables


def table_rows(number, caption, body):
    header, _, *rows = [
        [cell.strip().replace("*", "") for cell in line.strip("|").split("|")]
        for line in body.splitlines()
    ]
    if number == "1":
        header[0] = "Measure"
    assert all(len(row) == len(header) for row in rows), number
    if number == "S10":
        prefix = re.search(r"https://[^\s]+/(?=\.)", caption)[0]
        for row in rows:
            if row[0] == "10x Genomics":
                row[2] = prefix + row[2]
    return [header, *rows]


def section_test_rows(report):
    header = [
        "Chemistry", "Analysis", "Context", "Section or comparison",
        "Signal loss 32 µm", "Signal loss 64 µm", "Signal loss 256 µm",
        "Mito. elevation 32 µm", "Mito. elevation 64 µm", "Mito. elevation 256 µm",
        "100 UMIs/bin", "50 UMIs/bin", "25 UMIs/bin", "12 UMIs/bin", "Result",
    ]
    rows = [header]
    chemistry = analysis = context = ""
    # ponytail: Table S1 has several layouts. Keep paired-test summaries as text so
    # every count and qualifier survives. Export structured tests directly from 82
    # if individual summary values later need to be queried.
    for line in report.splitlines():
        value = line.strip()
        if not value or value.startswith("##") or line.startswith(("section  ", "test  ")):
            continue
        if line.startswith("# "):
            chemistry = value.removeprefix("# ").removesuffix(" sections")
            continue
        if line.startswith(("Neighbourhood test", "Section-wide test", "Region test",
                            "False-positive rate (%)")):
            analysis, context = value, ""
            continue
        if line.startswith("  ") and not line.startswith("   "):
            context = value.rstrip(":")
            continue

        cells = re.split(r" {2,}", value)
        best_bins = [""] * 6
        by_depth = [""] * 4
        if ": best bin size " in analysis and not line.startswith(" "):
            label, *best_bins = cells
            assert len(best_bins) == 6, line
            result = ""
        elif analysis.startswith("Region test") and " bins," in value:
            label, rest = value.split(":", 1)
            rates, result = rest.strip().split("; ", 1)
            by_depth = re.split(r" {2,}", rates)
            assert len(by_depth) == 4, line
        elif analysis.startswith("False-positive rate") and len(cells) == 7:
            test, size, *by_depth, result = cells
            label = f"{test}, {size}"
        else:
            label, result = re.split(r" {2,}|: ", value, maxsplit=1)
        assert chemistry and analysis, line
        rows.append([chemistry, analysis, context, label, *best_bins, *by_depth, result])
    assert len(rows) > 1 and all(len(row) == len(header) for row in rows)
    return rows


def export_final(tables, destination):
    table_csvs = {
        f"table{number}.csv": table_rows(number, caption, body)
        for number, (caption, body, _) in tables.items()
    }
    table_csvs["tableS1.csv"] = section_test_rows((OUT / "tableS1.txt").read_text(encoding="utf-8"))
    with TemporaryDirectory(dir=destination.parent) as directory:
        staging = Path(directory)
        for figure in FIGURES:
            source = ROOT / "methods_diagrams" if figure in DIAGRAMS else OUT
            name = f"{figure}.png"
            copyfile(source / name, staging / name)
        for name, rows in table_csvs.items():
            with (staging / name).open("w", encoding="utf-8", newline="") as file:
                csv.writer(file).writerows(rows)
        destination.mkdir(exist_ok=True)
        for path in staging.iterdir():
            path.replace(destination / path.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript", type=Path, default=ROOT / "docs" / "manuscript_draft.md",
                        help="main manuscript Markdown (default: docs/manuscript_draft.md)")
    parser.add_argument("--supplement", type=Path, default=ROOT / "docs" / "supplementary_draft.md",
                        help="supplementary Markdown (default: docs/supplementary_draft.md)")
    parser.add_argument("--final-output", action="store_true",
                        help="write only the final PNGs and CSVs to figures/final_output")
    args = parser.parse_args()
    for flag, draft in (("--manuscript", args.manuscript), ("--supplement", args.supplement)):
        if not draft.is_file():
            parser.error(f"{flag} requires an existing Markdown file. Supply it locally.")
    tables = manuscript_tables((args.manuscript, args.supplement))
    if args.final_output:
        export_final(tables, OUT / "final_output")
        return
    for number, (caption, body, note) in tables.items():
        (OUT / f"table{number}.md").write_text(f"{caption}\n\n{body}\n{note}", encoding="utf-8")


if __name__ == "__main__":
    main()
