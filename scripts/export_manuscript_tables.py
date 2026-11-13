#!/usr/bin/env python3
"""Export the manuscript tables and optionally package its final figures and tables."""
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


def export_final(tables, destination):
    table_csvs = {
        f"table{number}.csv": table_rows(number, caption, body)
        for number, (caption, body, _) in tables.items()
    }
    with TemporaryDirectory(dir=destination.parent) as directory:
        staging = Path(directory)
        for figure in FIGURES:
            source = ROOT / "methods_diagrams" if figure in DIAGRAMS else OUT
            name = f"{figure}.png"
            copyfile(source / name, staging / name)
        copyfile(OUT / "tableS1.txt", staging / "tableS1.txt")
        for name, rows in table_csvs.items():
            with (staging / name).open("w", encoding="utf-8-sig", newline="") as file:
                csv.writer(file).writerows(rows)
        destination.mkdir(exist_ok=True)
        for path in staging.iterdir():
            path.replace(destination / path.name)
        (destination / "tableS1.csv").unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript", type=Path, default=ROOT / "docs" / "manuscript_draft.md",
                        help="main manuscript Markdown (default: docs/manuscript_draft.md)")
    parser.add_argument("--supplement", type=Path, default=ROOT / "docs" / "supplementary_draft.md",
                        help="supplementary Markdown (default: docs/supplementary_draft.md)")
    parser.add_argument("--final-output", action="store_true",
                        help="write the final PNGs, CSVs and Table S1 TXT to figures/final_output")
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
