set -uo pipefail
WK=${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
. "$WK/scripts/00_slot.sh"
cd "$WK"
IMG=${IMG:-vhd-qc-benchmark}
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
py /work/scripts/69_figure_binsize.py --stack 100 50 25 12 --metric f1 --depth-titles --block-height 2.55 \
  --name fig4_binsize_f1 2>&1 | grep -c wrote | sed "s/^/--- 69 Fig. 4: files written /"
for s in scripts/63_figure_discarded.py scripts/91_figure_routes.py \
         scripts/95_figure_downstream.py scripts/90_supp_figs.py scripts/61_figures_grid.py; do
  echo "--- $s"
  py "/work/$s" 2>&1 | tail -3
done
py /work/scripts/compose_figures.py 2>&1 | tail -9

echo
echo "######## 4. generated tables"
py /work/scripts/93_k48_table.py > "$OUT/tableS3.md" 2>/dev/null &&
  echo "Table S3 written to $OUT/tableS3.md"
py /work/scripts/82_section_level_tests.py > "$OUT/table3_tests.txt" 2>&1 &&
  echo "Table 3 / S2 tests written to $OUT/table3_tests.txt"
grep -E '^#|sign test|rises in' "$OUT/table3_tests.txt"
py /work/scripts/126_manuscript_numbers.py 2>&1 | tail -1
python3 scripts/127_table1.py > "$OUT/table1.md" && echo "Table 1 written to $OUT/table1.md"
python3 scripts/137_table3.py 2>&1 | sed "s/^/Table 3 and S5: /"
py /work/scripts/129_boundary_corrections.py 2>&1 | tail -1
python3 scripts/128_guideline_figure.py 2>&1 | tail -2
python3 scripts/138_supp_numbers.py 2>&1 | tail -1
py /work/scripts/139_mito_threshold.py Visium_HD_Human_Kidney_FFPE_d100 Visium_HD_Mouse_Kidney_d100 2>&1 | tail -2
python3 scripts/117_methods_diagrams.py 2>&1 | tail -1

echo "### done $(date -Is)"
