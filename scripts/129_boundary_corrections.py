import importlib.util, os
import pandas as pd

import panel as P

GRID = "/work/benchmark_grid"
ARMS = ("09_mad_32um", "10_mad_64um", "34_mad_64um_percovered", "35_mad_64um_excludepartial",
        "41_mad_128um", "27_spotsweeper_64um_k48", "28_spotsweeper_128um_k48")
spec = importlib.util.spec_from_file_location("m79", "/work/scripts/79_build_model_tables.py")
m79 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m79)

rows = []
for sec, depth, path in P.cells(GRID):
    cell = os.path.basename(path)
    eb = m79.edge_barcodes(GRID, cell, 16.0)
    if eb is None:
        continue
    edge, n_tissue = eb
    for arm in ARMS:
        f = os.path.join(path, "results", arm, "clean.txt")
        if not os.path.exists(f):
            continue
        got = {ln.strip() for ln in open(f) if ln.strip()}
        rows.append(dict(section=sec, chem=P.chemistry(sec), depth=depth, arm=arm,
                         edge_disc=len(got & edge), edge_n=len(edge),
                         disc=len(got), n_tissue=int(n_tissue)))
out = pd.DataFrame(rows)
p = os.path.join(GRID, "figures", "boundary_corrections.csv")
out.to_csv(p, index=False)
print("wrote", p, len(out), "rows")
