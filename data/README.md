# data

| File | Content |
|---|---|
| `sources.tsv` | The 40 sections: species, chemistry, provider, download URL, licence and the study to cite |
| `big512_runs.txt` | The 67 runs that `scripts/run_section.sh` gives 512 µm LQRs and 128 µm bins |
| `widened_runs.txt` | The 8 runs on which the widened neighbourhoods (192, 768 and 3072 neighbours) are run |
| `reconstructed/` | The 28 injector records rebuilt for the 24 runs whose first records were overwritten (below) |

## How the deposited runs were built

`scripts/run_section.sh` builds a run in one pass, as the Methods describe. The deposited runs
were built by the same programs in stages, as the design grew. Each stage is one call of the
injector, `scripts/10_inject_artifacts.py`, and each call leaves a record of its arguments and
seed (`injection_params*.json`) and a manifest of the conditions it wrote (`manifest*.csv`).
Three things about the stages matter for an exact rebuild.

**The mitochondrial rule.** The 109 runs built before 2026-10-14 first carried mitochondrial
LQRs on an absolute target of 25% of each bin's counts. Their mitochondrial conditions were
injected again on the relative rule (each section's baseline odds raised tenfold) by
`scripts/118_reinject_mito.py`, which calls the injector and leaves its own records,
`injection_params_mito_relative_{core,edge,big}.json`. The injector draws every condition of a
call from one random stream, in the order the conditions are written, so in these runs the
signal-loss conditions of 64 and 256 µm follow the draws of the old mitochondrial conditions.
They are reproduced by replaying the first call, with its absolute target, and the
mitochondrial conditions by replaying the 118 call.

**Edge and 512 µm conditions.** On 17 runs of four pilot sections (listed below), the first
edge conditions were added from 2026-09-21 to 09-23, and on human pancreas and mouse embryo at
100 and 12 UMIs per bin the first 512 µm conditions on 2026-09-22, by calls without a record
suffix. Every other run received its edge conditions from 2026-10-15 to 10-23, and the runs in
`big512_runs.txt` their 512 µm LQRs on 2026-10-25 and 26, each by a call with its own record
(`_edge`, `_big512`). The four early 512 µm runs are judged by every arm.

**Overwritten records.** On 24 runs of the pilot sections (23 of the paper's 155 runs and mouse
embryo at 200 UMIs per bin), a later call into the run was made without a record suffix and
wrote over `injection_params.json`, the record of the call that built the run: a trial of a
second mitochondrial rule on 7 runs (2026-09-20), and the first edge or 512 µm call on 17
(2026-09-21 to 09-23). On human pancreas at 100 and 12 UMIs per bin the 512 µm call also wrote
over the record of the edge call, and on mouse embryo at 100 and 12 the edge call over that of
the 512 µm call. The conditions those calls wrote are still the ones on disk, apart from the
mitochondrial conditions that 118 replaced. `reconstructed/` gives each of the 28 lost records,
rebuilt from the command that made it; its `reconstructed` field says which.

| Runs | `injection_params.json` now holds | Reconstructed |
|---|---|---|
| Human kidney FFPE at 100, 50, 25 and 12; mouse kidney at 100; mouse small intestine at 100 and 12 UMIs per bin | the trial of a second mitochondrial rule (kind label `mito=mitorel`), whose conditions are no longer on disk | the build call |
| Human colon cancer and mouse brain at 100, 50, 25 and 12; human pancreas at 50 and 25; mouse embryo at 200, 50 and 25 | the first edge call | the build call |
| Mouse embryo at 100 and 12 | the first edge call | the build call and the first 512 µm call (`_big512`) |
| Human pancreas at 100 and 12 | the first 512 µm call | the build call and the first edge call (`_edge`) |

## Verification

`scripts/136_rebuild_check.py` replays records on the public section and compares each
condition they write with the deposited one: the 8 µm matrix, its barcodes, the bin positions
and the ground truth. With `ONLY_GIVEN=1` it replays only the records given and checks only the
conditions they write.

Each reconstructed record was replayed with `ONLY_GIVEN=1`, together with the run's 118 record
of the mitochondrial conditions it shares, on at least one run per kind of call:

| Reconstructed call | Runs replayed | Result |
|---|---|---|
| the build call | human colon cancer and human pancreas at 100, 50, 25 and 12; human kidney FFPE at 50, 25 and 12 UMIs per bin | the undamaged and the three signal-loss conditions reproduced exactly on all 11 runs |
| the first edge call | human pancreas at 12 | the three edge signal-loss conditions reproduced exactly |
| the first 512 µm call | mouse embryo at 12 | the 512 µm signal-loss condition reproduced exactly |

The mitochondrial conditions these calls wrote under the absolute rule were later replaced, and
most of them differ from the deposited ones, as they should; the 118 records replayed alongside
reproduced every deposited mitochondrial condition they cover exactly. The other 15
reconstructed records differ from the replayed ones only in depth or section and were not
replayed; their `reconstructed` field says so.

On a run whose matrices are on disk, replay a reconstructed record with the 118 record that
replaced its mitochondrial conditions:

    ONLY_GIVEN=1 bash scripts/rebuild_run.sh Visium_HD_Human_Kidney_FFPE_d100 \
      /work/data/reconstructed/Visium_HD_Human_Kidney_FFPE_d100.json \
      /work/benchmark_grid/Visium_HD_Human_Kidney_FFPE_d100/injection_params_mito_relative_core.json
