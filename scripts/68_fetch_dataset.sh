set -euo pipefail
NAME=$1
VER=${2:-4.0.1}
WK=${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
DEST=$WK/public/$NAME/binned_outputs
CF=https://cf.10xgenomics.com/samples/spatial-exp
S3=https://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp
URL=$CF/$VER/$NAME/${NAME}_binned_outputs.tar.gz
curl -sfI "$URL" >/dev/null || URL=$S3/$VER/$NAME/${NAME}_binned_outputs.tar.gz
[ -d "$DEST/square_008um/raw_feature_bc_matrix" ] && { echo "$NAME already present"; exit 0; }
TMP=$(mktemp -d "$WK/public/.fetch_${NAME}_XXXX"); trap 'rm -rf "$TMP"' EXIT
echo "$NAME: streaming $URL"
curl -sSL --retry 3 --retry-delay 10 "$URL" \
  | tar -xz -C "$TMP" --wildcards '*/square_008um/filtered_feature_bc_matrix/*' \
                      '*/square_008um/raw_feature_bc_matrix/*' \
                      '*/square_008um/spatial/tissue_positions.parquet' \
                      '*/square_008um/spatial/scalefactors_json.json'
SRC=$(dirname "$(find "$TMP" -type d -name square_008um -print -quit)")
[ -n "$SRC" ] || { echo "$NAME: square_008um not found in tarball"; exit 1; }
for f in filtered_feature_bc_matrix/matrix.mtx.gz raw_feature_bc_matrix/matrix.mtx.gz \
         spatial/tissue_positions.parquet spatial/scalefactors_json.json; do
  [ -s "$SRC/square_008um/$f" ] || { echo "$NAME: missing $f"; exit 1; }
done
n=$(zcat "$SRC/square_008um/raw_feature_bc_matrix/barcodes.tsv.gz" | wc -l)
k=$(zcat "$SRC/square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz" | wc -l)
[ "$n" -gt "$k" ] || { echo "$NAME: raw $n <= filtered $k"; exit 1; }
mkdir -p "$DEST" && mv "$SRC/square_008um" "$DEST/square_008um"
echo "$NAME: in place — $k in-tissue of $n capture-area bins"
