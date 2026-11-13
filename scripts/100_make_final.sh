#!/usr/bin/env bash
set -uo pipefail
WK=${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
. "$WK/scripts/00_slot.sh"
cd "$WK"
IMG=${IMG:-vhd-qc-benchmark:0.4}
OUT=benchmark_grid/figures
py() { wait_for_slot >&2; docker run --rm -v "$WK":/work -w /work "$IMG" /opt/pyenv/bin/python "$@"; }

echo "######## 1. completeness gate"
LEFT=$(python3 scripts/104_leftovers.py | grep -c .)
echo "cells still needing work: $LEFT (must be 0)"
[ "$LEFT" -eq 0 ] || python3 scripts/104_leftovers.py | cut -c1-120
echo

echo "######## 2. inputs the figures read"
py /work/scripts/79_build_model_tables.py --outdir /work 2>&1 | tail -3
py /work/scripts/85_region_ci.py 2>&1 | tail -1

echo "######## 3. figures"
py /work/scripts/86_figure_tool.py 2>&1 | tail -3
for d in 100 50 25 12; do for m in recall f1; do
  py /work/scripts/69_figure_binsize.py --depth "$d" --metric "$m" 2>&1 | grep -c wrote |
    sed "s/^/--- 69 d$d $m: files written /"
done; done
py /work/scripts/69_figure_binsize.py --stack 100 12 --metric recall --name fig2_binsize_recall 2>&1 |
  grep -c wrote | sed "s/^/--- 69 Fig. 2: files written /"
py /work/scripts/69_figure_binsize.py --stack 100 50 25 12 --metric f1 --overlay 256 \
  --name fig5_binsize_f1 2>&1 | grep -c wrote | sed "s/^/--- 69 Fig. 5: files written /"
for s in scripts/63_figure_discarded.py scripts/91_figure_routes.py \
         scripts/95_figure_downstream.py scripts/90_supp_figs.py scripts/61_figures_grid.py; do
  echo "--- $s"
  py "/work/$s" 2>&1 | tail -3
done
demo_sections=(Visium_HD_Human_Colon_Cancer Visium_HD_Human_Kidney_FFPE
               Visium_HD_Human_Lung_Cancer Visium_HD_3prime_Human_Ovarian_Cancer)
for section in "${demo_sections[@]}"; do
  if [[ ! -s "$OUT/demo_downstream/$section.png" ]]; then
    py /work/scripts/122_demo_downstream_maps.py --sections "$(IFS=,; printf '%s' "${demo_sections[*]}")" || exit 1
    break
  fi
done
py /work/scripts/146_molecule_curve.py > logs/146_molecule_curve.log || exit 1
py /work/scripts/compose_figures.py 2>&1 | tail -9

echo
echo "######## 4. generated tables"
py /work/scripts/93_k48_table.py > "$OUT/tableS7_neighbours_generated.md" 2>/dev/null &&
  echo "Table S7 source written to $OUT/tableS7_neighbours_generated.md"
py /work/scripts/149_sa_cutoff.py 2>&1 | grep -E 'bins, 12 UMIs'
py /work/scripts/82_section_level_tests.py > logs/82_section_level_tests.log 2>&1 &&
  echo "Table S1 written to $OUT/tableS1.txt"
grep -E '^#|sign test|rises in' "$OUT/tableS1.txt"
py /work/scripts/126_manuscript_numbers.py 2>&1 | tail -1
python3 scripts/127_table1.py > "$OUT/tableS8_panel_generated.md" && echo "Table S8 source written to $OUT/tableS8_panel_generated.md"
python3 scripts/137_table3.py 2>&1 | sed "s/^/Table S2 source: /"
py /work/scripts/129_boundary_corrections.py 2>&1 | tail -1
python3 scripts/128_guideline_figure.py 2>&1 | tail -2
python3 scripts/138_supp_numbers.py 2>&1 | tail -1
py /work/scripts/139_mito_threshold.py Visium_HD_Human_Kidney_FFPE_d100 Visium_HD_Mouse_Kidney_d100 2>&1 | tail -2
python3 scripts/117_methods_diagrams.py 2>&1 | tail -1

echo
echo "######## 5. numbers still written by hand in the draft"
echo "These were measured at the old neighbour count and must be checked against what"
echo "section 3 just regenerated. See the latest docs/progress_*.md for what each is."
grep -nP '\b(2304|576|144)\b(?!,\d{3})|= 36\b' docs/manuscript_draft.md | cut -c1-100
echo
echo "######## 6. invariants"
printf 'em-dashes:      %s main, %s supplementary (must be 0)\n' \
  "$(grep -c '—' docs/manuscript_draft.md)" "$(grep -c '—' docs/supplementary_draft.md)"
printf 'citation holes: %s (must be 0)\n' "$(grep -c '<citation>' docs/manuscript_draft.md)"
printf 'open TODO markers: %s main, %s supplementary (must be 0 at submission)\n' \
  "$(grep -c 'TODO' docs/manuscript_draft.md)" "$(grep -c 'TODO' docs/supplementary_draft.md)"
printf 'words:          %s main, %s supplementary\n' \
  "$(wc -w < docs/manuscript_draft.md)" "$(wc -w < docs/supplementary_draft.md)"
echo "### done $(date -Is)"
