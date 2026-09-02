# PCA arm crossover: CPU vs GPU

Wall-clock and phase timings for the PCA stage across four arms, three dataset
sizes and two machines. `benchmark_pca_crossover.yaml` at the repo root is the
plan that produced them.

## Running it

    ob run benchmark_pca_crossover.yaml --cores 8 --until PCA -k \
       --with-capability gpu --unpinned --out-dir out_main   # GPU box
    ob run benchmark_pca_crossover.yaml --cores 8 --until PCA -k \
       --out-dir out_main                                    # mac / any CPU box

`pc-rapids` declares `requires_capabilities: [gpu]`, so without
`--with-capability gpu` it is pruned automatically -- there is no separate
CPU-only variant of this plan to keep in sync.

`--unpinned` is only needed for the GPU run: `pc-rapids` pins the branch
`testing`, so its resolved SHA is whatever HEAD was at run time.

### omnibenchmark version is load-bearing

`pixi.toml` pins omnibenchmark to `1fddceb` (0.6.0). Three things this
benchmark depends on landed after 0.5.3:

* `6063db9` emits snakemake's `threads:` from `resources.cores`. Without it the
  PCA stage's `resources: cores: 8` is inert: every job is scheduled as
  1 thread, so `--cores 8` runs 7 PCA arms concurrently and they contend.
  Measured on tenx-0020k, that inflated pc-scanpy from 5.4s to 68.3s.
* `--until STAGE`, which prunes before module resolution -- so CNTFCT
  (pinned to a branch) does not have to resolve for a PCA-only run.
* capability gating, which is what lets one plan serve both machines.

### The tenx rung is machine-local

The `tenx-0020k` DATA module points at a `file://` path on the GPU box. On
another machine, repoint that one `uri:` field. be1 and d-pca are fetched, so
they need no edit.

## Files

* `performances_gpu.tsv` -- 162 records, Ryzen 9 7940HS + RTX 2000 Ada
* `performances_mac.tsv` -- 117 records, Apple M4 (10 cores)
* `phases_gpu.tsv` -- 315 phase records (load/compute/write) from the two
  instrumented modules

## What the numbers say

Medians of 5 seeds, tenx-0020k (19,696 cells):

| arm | GPU box | mac M4 |
|-----|---------|--------|
| pc-rapids halko  | 5.40 | (gated out) |
| pc-scanpy arpack | 6.59 | 5.44 |
| pc-rcppml lanczos| 7.29 | 4.16 |
| pc-scrapper irlba| 11.47| 6.75 |

Wall clock understates the GPU badly. Splitting by phase at 19,696 cells:

| | load | compute | write | import/startup | wall |
|---|---|---|---|---|---|
| pc-rapids | 0.52 | **0.77** | 0.94 | 3.02 | 5.40 |
| pc-scanpy | 0.51 | **4.13** | 0.07 | 1.88 | 6.59 |

On compute alone rapids is 5.35x faster, and the ratio grows with n
(0.90x at 1.7k, 2.65x at 3.8k, 5.35x at 19.7k). The compute-only crossover is
around 1.8k cells; the wall-clock crossover is around 9.9k. The gap between
those two numbers is harness overhead, and it is where the tuning is:

1. `pc-rapids` writes its TSV with pandas `to_csv` (0.94s) where `pc-scanpy`
   uses polars `write_csv` (0.07s) for the same output -- 13x, and more than
   rapids' entire GPU compute.
2. Process start + imports is 1.9-3.0s per job, larger than load+compute+write
   combined for rapids. This is the dominant fixed cost.
3. `load` is 27% of rapids' instrumented time and grows linearly; an mmap'd
   binary input would cut it.

## Caveats

* `max_rss` is recorded on Linux (449-1987 MB across arms) but comes back 0 on
  macOS, so the RAM axis is Linux-only.
* Only `pc-rapids` and `pc-scanpy` emit phases. The R arms (`pc-scrapper`,
  `pc-rcppml`) have no instrumentation, and `pc-rcppml` on the M4 is the one
  arm rapids does not beat -- so the most interesting comparison is the one
  that cannot be decomposed. `phases.py` in the rapids module says it should be
  promoted to obkit; doing that and calling it from the R modules would close
  this.
* Three dataset sizes spanning one decade, all dominated by fixed costs. The
  slopes here are a local tangent, not a scaling law.
* Cores were fixed at 8 everywhere. That axis is unsampled.
