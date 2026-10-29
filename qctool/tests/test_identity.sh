set -uo pipefail
T=${1:-Visium_HD_Mouse_Embryo_d12}; C=${2:-dropout01edge_64um}
B=/work/benchmark_grid/$T; Q=/work/qctool; TMP=$(mktemp -d)
fail=0
check() {
  if [ ! -f "$3" ]; then echo "SKIP $1: no $3"; return; fi
  if diff -q <(sort -u "$2") <(sort -u "$3") > /dev/null; then echo "PASS $1: $(sort -u "$2" | wc -l) bins, identical to the arm"
  else echo "FAIL $1: tool $(sort -u "$2" | wc -l), arm $(sort -u "$3" | wc -l), differ by $(comm -3 <(sort -u "$2") <(sort -u "$3") | wc -l)"; fail=1; fi
}
while read -r name test bin k arm; do
  s=$(date +%s)
  Rscript $Q/qc_run.R --views $B/$C --test $test --bin $bin --k $k --source8 /src/square_008um \
    --out $TMP/$name.txt > $TMP/$name.log 2>&1 || { echo "FAIL $name: qc_run.R error"; tail -3 $TMP/$name.log; fail=1; continue; }
  check "$name ($(( $(date +%s) - s )) s)" $TMP/$name.txt $B/results/$arm/$C.txt
done <<'L'
mad_8 mad 8 48 02_mad_8um
mad_32 mad 32 48 09_mad_32um
ss_16_k48 spotsweeper 16 48 25_spotsweeper_16um_k48
ss_64_k48 spotsweeper 64 48 27_spotsweeper_64um_k48
sa_8 spatialartifacts 8 48 19_spatialartifacts_8um
sa_16 spatialartifacts 16 48 20_spatialartifacts_16um
ss_8_k192 spotsweeper 8 192 44_spotsweeper_8um_k192
L

mkdir -p $TMP/in/square_008um/spatial
ln -s $B/$C/square_008um/filtered_feature_bc_matrix $TMP/in/square_008um/
cp /src/square_008um/spatial/tissue_positions.parquet /src/square_008um/spatial/scalefactors_json.json $TMP/in/square_008um/spatial/
/opt/pyenv/bin/python - "$TMP" "$B/$C" <<'PY' || fail=1
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, "/work/qctool"); import vhdqc
tmp, cond = sys.argv[1], sys.argv[2]
mat, bc8, feat = vhdqc.inj.read_mex(f"{tmp}/in/square_008um/filtered_feature_bc_matrix")
vhdqc.prepare(f"{tmp}/in", f"{tmp}/views", [16, 32, 64, 128], mat, bc8, feat)
bad = 0
for b in (8, 16, 32, 64, 128):
    v, r = f"{tmp}/views/square_{b:03d}um", f"{cond}/square_{b:03d}um"
    if not os.path.exists(r):
        print(f"SKIP view {b} um: not in the condition"); continue
    m1, c1, _ = vhdqc.inj.read_mex(f"{v}/filtered_feature_bc_matrix")
    m2, c2, _ = vhdqc.inj.read_mex(f"{r}/filtered_feature_bc_matrix")
    p1 = pd.read_parquet(f"{v}/spatial/tissue_positions.parquet")
    p2 = pd.read_parquet(f"{r}/spatial/tissue_positions.parquet")
    ok = (list(c1) == list(c2) and m1.shape == m2.shape and (m1 != m2).nnz == 0
          and all((p1[k].to_numpy() == p2[k].to_numpy()).all() for k in ("barcode", "array_row", "array_col",
                                                                          "n_children", "in_tissue", "is_truncated"))
          and np.allclose(p1[["pxl_row_in_fullres", "pxl_col_in_fullres"]], p2[["pxl_row_in_fullres", "pxl_col_in_fullres"]], atol=1e-6))
    print(("PASS" if ok else "FAIL") + f" view {b} um: {len(c1):,} bins")
    bad |= not ok
sys.exit(bad)
PY
rm -rf "$TMP"
[ $fail -eq 0 ] && echo "ALL PASS" || echo "SOME FAILED"
