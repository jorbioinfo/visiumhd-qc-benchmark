# visiumhd-qc-benchmark

Code for the paper *How should sub-cellular spatial transcriptomics be quality controlled? A guideline from a benchmark with injected damage*. It builds a benchmark with known damage on 40 public Visium HD sections,
runs each quality control (QC) method on it, scores every method against the damage, and
regenerates every figure, table and number in the paper. It also holds the implementation
of aggregation and deconvolution (A&D) and a command-line tool that applies the paper's QC
procedure to a new section.

The benchmark data (injection records, ground truth, every arm's verdict and the scored
tables) are in Zenodo: https://doi.org/10.5281/zenodo.23681225.

## What the benchmark does

1. **Depth control.** Each section is thinned by binomial sampling to a mean of 100, 50, 25
   and 12 UMIs per 8 µm bin. A run is one section at one depth.
2. **Synthetic low-quality regions (LQRs).** Circles of 32, 64 and 256 µm are injected at
   known positions as signal loss (each count kept with probability 0.01) or mitochondrial
   elevation (each section's baseline odds raised tenfold), inside the tissue and at its
   edge. On the sections large enough to hold them, 512 µm LQRs are added.
3. **QC arms.** Each arm is one QC method at one bin size: MAD thresholds (scuttle),
   SpotSweeper, SpatialArtifacts, miQC, a fixed UMI cutoff and no filter. An arm on bins
   coarser than 8 µm runs on summed bins and projects each verdict back to the 8 µm bins
   inside (A&D), so every arm reports on the same 8 µm bins.
4. **Scoring.** Recall, cost on the uninjected condition, F1, a marker check of what was
   removed, and the agreement of downstream clusterings between depths.

## Layout

| Path | Content |
|---|---|
| `docker/Dockerfile` | The environment every step runs in |
| `R/ad.R` | A&D: the projection of coarse verdicts back to 8 µm bins |
| `scripts/` | Every step of the benchmark, numbered as they were written |
| `scripts/run_section.sh` | Build and score one run as the Methods describe |
| `scripts/rebuild_run.sh` | Replay one run's injection records and compare with its matrices |
| `scripts/100_make_final.sh` | Regenerate every figure, table and number from a complete grid |
| `qctool/` | The command-line QC tool (the procedure of Fig. 6) |
| `data/` | The panel's sources, which runs carry the subset conditions, and the reconstructed injection records (see `data/README.md`) |
| `checks/` | Minimal examples of two package issues the benchmark works around |

The main steps, in the order a run passes them:

| Script | Step |
|---|---|
| `68_fetch_dataset.sh`, `113_zenodo_extract.py` | Download a section's 8 µm files (10x Genomics, Zenodo) |
| `10_inject_artifacts.py` | Thin to a depth, place the LQRs, write the conditions and their coarse views |
| `75_add_coarse_bin.py` | Add 128 µm views to a run |
| `12_run_arms_k48all.R` | Run the QC arms |
| `11_score_benchmark.py` | Score every arm on every condition |
| `50_marker_validation.py` | Marker check of the removed bins |
| `94_downstream_binned.py` | Clustering agreement between 100 and 12 UMIs per bin |
| `136_rebuild_check.py` | Replay a run's recorded injector calls and compare with the run |

The figure and table scripts that `100_make_final.sh` calls are listed in it, each with what
it reads and writes in its header.

## Requirements

Docker, about 60 GB of memory for the largest sections (the 11 mm slide), and disk for the
runs: 1.5 to 9 GB each, about 630 GB for the 155 runs of the paper. Build the image from the repository root:

    docker build -f docker/Dockerfile -t vhd-qc-benchmark .

The benchmark ran in an image built from this file on 2026-09-23 on Bioconductor 3.23:
R 4.6.1, scuttle 1.22.0, SpotSweeper 1.8.0, SpatialArtifacts 1.0.0, miQC 1.20.0,
SpatialExperiment 1.22.0; Python 3.12.3 with NumPy 2.1.3, SciPy 1.14.1, pandas 2.2.3,
pyarrow 18.1.0, scanpy 1.11.5, anndata 0.11.4, leidenalg 0.10.2, igraph 0.11.9 and
scikit-learn 1.5.2. All packages are release builds used as published.

Every container goes through `scripts/00_slot.sh`, which holds the number of benchmark
containers at `CAP` (default 5). Paths are relative to the repository root, or to
`VHDQC_ROOT` when it is set; inside the container the root is mounted at `/work` and the
section's `binned_outputs` at `/src`.

## Running the benchmark

Fetch a section, then build its runs:

    bash scripts/68_fetch_dataset.sh Visium_HD_Mouse_Embryo 4.0.1     # name and Space Ranger version, see data/sources.tsv
    for d in 100 50 25 12; do bash scripts/run_section.sh Visium_HD_Mouse_Embryo $d; done

The six-cancer sections are read from their Zenodo zip without downloading it whole:

    python3 scripts/113_zenodo_extract.py visiumhd_Bladder.zip Visium_HD_SixCancer_Bladder_Cancer 297675

Then the downstream clustering, one section at a time (100 against 12 UMIs per bin):

    docker run --rm -v "$PWD":/work vhd-qc-benchmark /opt/pyenv/bin/python \
      /work/scripts/94_downstream_binned.py --section Visium_HD_Mouse_Embryo --bin 32 \
      --arms 01_no_qc,02_mad_8um,10_mad_64um,03_spotsweeper_8um,24_spotsweeper_8um_k48

and, once every run is scored, every figure, table and number:

    bash scripts/100_make_final.sh

## Reproducing the deposited runs exactly

The deposit holds each run's injection records (`injection_params*.json`), one per call of
the injector, and the ground truth of every condition, but not the count matrices. A call is
deterministic given its arguments and seed, so replaying a run's calls on the public section,
in the order the run was built (the call that built it, then the edge, 512 µm and
mitochondrial re-injection calls), rebuilds its matrices bin for bin. To rebuild a run, call
`scripts/10_inject_artifacts.py` with each record's arguments in that order;
`scripts/136_rebuild_check.py` shows how a record maps to the arguments. Where the run's
matrices are on disk, `scripts/rebuild_run.sh` replays its calls and compares every condition
with them:

    bash scripts/rebuild_run.sh Visium_HD_Human_Lymph_Node_FFPE_d12

The deposited runs were built in stages as the design grew, so not every run was made in
the single pass of `run_section.sh`. On 24 runs the record of the call that built the run was
overwritten by a later call. `data/README.md` describes the stages, lists these runs and gives
their reconstructed records, in `data/reconstructed/`, with how they were checked.

## The QC tool

`qctool/vhdqc.py` applies the procedure of Fig. 6 to a Space Ranger output. It asks for the
two things the data cannot give, the size of any damage seen on the tissue image and
whether it sits at the tissue edge, and reads the test, bin size and number of neighbours
from `qctool/guideline_settings.json`, the file the script that draws Fig. 6 writes. Each
test is run by the function the benchmark arm ran (`qctool/qc_core.R`, cut verbatim from
`scripts/12_run_arms_k48all.R` and `R/ad.R` by `qctool/make_core.R`).

    docker run --rm -v "$PWD":/work -v /path/to/outs/binned_outputs:/in:ro vhd-qc-benchmark \
      /opt/pyenv/bin/python /work/qctool/vhdqc.py --input /in --out /work/qc_out --damage-um 200

`qctool/tests/test_identity.sh` checks that the tool's verdicts equal the benchmark arms'
on a deposited run.

## Citation

[citation of the paper, to be added on publication]

## Licence

MIT, see `LICENSE`. The Visium HD sections are released by their producers under CC BY 4.0
(`data/sources.tsv`).
