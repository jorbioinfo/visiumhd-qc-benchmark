import gzip, io, json, os, shutil, sys, tempfile, time, urllib.request, zipfile

def fetch(req, timeout):
    for attempt in range(6):
        try:
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception as e:
            if attempt == 5: raise
            wait = 30 * 2 ** attempt
            print(f"  retry in {wait} s after: {e}", flush=True)
            time.sleep(wait)

KEY, NAME, BINS = sys.argv[1], sys.argv[2], int(sys.argv[3])
WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = f"{WK}/public/{NAME}/binned_outputs/square_008um"
if os.path.isdir(f"{DEST}/raw_feature_bc_matrix"):
    sys.exit(print(f"{NAME} already present") or 0)

class RangeFile(io.RawIOBase):
    def __init__(s, url, size): s.url, s.size, s.pos = url, size, 0
    def seekable(s): return True
    def readable(s): return True
    def tell(s): return s.pos
    def seek(s, off, wh=0):
        s.pos = off if wh == 0 else s.pos + off if wh == 1 else s.size + off; return s.pos
    def readinto(s, b):
        if s.pos >= s.size: return 0
        end = min(s.pos + len(b), s.size) - 1
        d = fetch(urllib.request.Request(s.url, headers={"Range": f"bytes={s.pos}-{end}"}), 600)
        b[:len(d)] = d; s.pos += len(d); return len(d)

rec = json.loads(fetch("https://zenodo.org/api/records/17999961", 120))
f = next(x for x in rec["files"] if x["key"] == KEY)
z = zipfile.ZipFile(io.BufferedReader(RangeFile(f["links"]["self"], f["size"]), buffer_size=64 << 20))
need = {}
for n in z.namelist():
    tail = n.split("binned_outputs/square_008um/")[-1] if "binned_outputs/square_008um/" in n else None
    if tail in ("spatial/tissue_positions.parquet", "spatial/scalefactors_json.json") or (
            tail and tail.split("/")[0] in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix")
            and tail.endswith((".mtx.gz", ".tsv.gz"))):
        need[tail] = n
want = [f"{m}/{x}" for m in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix")
        for x in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")] + \
       ["spatial/tissue_positions.parquet", "spatial/scalefactors_json.json"]
missing = [w for w in want if w not in need]
if missing:
    sys.exit(f"{NAME}: missing in {KEY}: {missing}")
tmp = tempfile.mkdtemp(prefix=f".fetch_{NAME}_", dir=f"{WK}/public")
try:
    for tail in want:
        out = os.path.join(tmp, "square_008um", tail)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with z.open(need[tail]) as src, open(out, "wb") as dst:
            shutil.copyfileobj(src, dst, 16 << 20)
        print(f"  {tail}  {os.path.getsize(out) / 1e6:.1f} MB", flush=True)
    k = sum(1 for _ in gzip.open(os.path.join(tmp, "square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz"), "rt"))
    n = sum(1 for _ in gzip.open(os.path.join(tmp, "square_008um/raw_feature_bc_matrix/barcodes.tsv.gz"), "rt"))
    if k != BINS:
        sys.exit(f"{NAME}: MISMATCH {k} in-tissue bins, web_summary reports {BINS}")
    if n <= k:
        sys.exit(f"{NAME}: raw {n} <= filtered {k}")
    os.makedirs(os.path.dirname(DEST), exist_ok=True)
    shutil.move(os.path.join(tmp, "square_008um"), DEST)
    print(f"{NAME}: in place, {k} in-tissue of {n} capture-area bins")
finally:
    shutil.rmtree(tmp, ignore_errors=True)
