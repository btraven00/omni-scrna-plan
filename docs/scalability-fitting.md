# Fitting scalability on the TENx rungs

Everything learned the hard way while getting a first cost ladder to run. Read
with [`nondominated-pipelines.md`](nondominated-pipelines.md), which says how a
slope enters the dominance framework (short version: it does **not** get an
epsilon — it sets the cost coordinate, and only `cost(n*) = a·n*^b` goes on the
front).

The slice is [`slices/scalability.yaml`](../slices/scalability.yaml); the fits
are `omni-scrna-analysis/fit_scaling.py`, the plot `plot_scaling.py`.

## 1. The model, and what each half is for

    cost(n) = a · n^b        log cost = log a + b log n

* **b (slope)** should be characteristic of the algorithm and travel between
  machines.
* **a (intercept)** is machine-, implementation- and language-dependent.
  Measured here: the two Leiden-on-a-graph arms have intercepts 6x apart
  (2.09e-6 vs 3.48e-7) with slopes within 5% (1.251 vs 1.316).

Only the **pair** predicts a cost, so report both plus the fit quality, and
state the n\* any prediction is evaluated at.

## 2. R² is not the diagnostic. The slope CI is.

With 4 rungs there are 2 residual df and t(0.975, 2) = 4.30. Measured:

| arm | slope | 95% CI | R² |
|---|---|---|---|
| CLUST/cl-scanpy | 1.251 | [1.03, 1.47] | 0.9967 |
| PCA/pc-scrapper | 0.974 | [0.73, 1.22] | 0.9933 |
| PCA/pc-scanpy | 0.879 | [0.59, 1.17] | 0.9885 |
| CLUST/cl-scrapper | 1.316 | [0.83, 1.80] | 0.9854 |
| CLUST/cl-seurat | 0.841 | [0.37, 1.31] | 0.9669 |
| NNG/nn-scanpy | 0.544 | [-0.65, 1.74] | 0.6566 |

R² above 0.98 with a CI spanning 0.59–1.17 is not a good fit; it is four points
on a monotone curve, where high R² is close to automatic. **Never quote a slope
without its CI.** `cl-seurat`'s interval spans sublinear to superlinear, so the
tempting claim "seurat scales best and overtakes scanpy at ~1e7 cells" is not
supportable from this data. It was made, and it was wrong.

## 3. Repeats per rung buy more than extra rungs

df = r·k − 2 for r repeats at each of k rungs.

| design | points | df | t(0.975) |
|---|---|---|---|
| 4 rungs × 1 | 4 | 2 | 4.30 |
| 6 rungs × 1 | 6 | 4 | 2.78 |
| 6 rungs × 3 | 18 | 16 | 2.12 |

Repeats at the small rungs are nearly free and are the only way to separate
measurement noise from a real kink. Do them before reaching for bigger rungs.

## 4. Measurement discipline

* **Serial, on an idle machine.** `resources: cores: N` + `--cores N` makes
  snakemake reserve the whole box per rule, so one job runs at a time. The
  accuracy sweeps ran `--cores 4` with four jobs in flight; those timings are
  not fit for a cost claim and must not be mixed in.
* **Cold cache.** `io_in` currently reads 0.0 on every node — that is the page
  cache, not the absence of reads. Drop caches between rungs or the I/O slopes
  are fiction.
* **Stable clocks.** Fixed CPU governor, turbo off. Drift across a multi-hour
  run reads as slope.
* **One resource, one fit.** Time and RAM have different exponents; fit
  separately. Disk too, once an out-of-core arm is in play.

## 5. Nothing is running on 8 cores (verify, don't assume)

`resources: cores: 8` reached snakemake (the generated Snakefile carries
`threads:` on 99 rules) and the algorithms still did not use it:

| arm | cpu_time/wall | mean_load% | effective cores |
|---|---|---|---|
| CLUST/cl-scanpy | 0.87 | 71 | ~0.9 |
| CLUST/cl-seurat | 0.94 | 94 | ~0.9 |
| NNG/nn-scanpy | 0.81 | 77 | ~0.8 |
| PCA/pc-scanpy | 1.21 | 100 | ~1.2 |
| PCA/pc-scrapper | 1.21 | 121 | ~1.2 |

Leiden on a graph is serial; sparse ARPACK is a sequential Lanczos whose matvec
is not a BLAS-3 path (scanpy's own `pca.py` header says so); `--num_threads 8`
into libscran's irlba buys 0.2 of a core. So `cores: 8` bought clean serial
execution and wasted 7 cores.

**Check this every time** with `cpu_time / s` and `mean_load` from
`performance.txt`. A sublinear slope was once explained here as 8-thread Amdahl
efficiency; that explanation was wrong because the stage was single-threaded,
and the real cause of `nn-scanpy`'s 0.544 with R² 0.66 is still unidentified —
pynndescent's index-build vs query balance shifting with n is the candidate.
It is also the most expensive stage (~20s at 39k, 4x anything else), so it is
the one most worth explaining.

## 6. Every arm must load the same way

`pc-scrapper` used to open the matrix with `TENxMatrix` (a DelayedArray) and
then immediately coerce to `dgCMatrix`: it paid block-wise HDF5 read cost
(1.77s vs scanpy's 0.16s at 39k, 11x) **and materialised anyway**, so it got
none of the out-of-core memory benefit. Its peak-RSS slope was the steepest of
any arm (0.255) — the opposite of what streaming looks like.

Now explicit: `--backend memory` (direct CSC read via rhdf5, matched to scanpy
and rapids) or `--backend delayed` (the old path, for the top rungs where
materialising is impossible). **`delayed` still coerces**, so it is a read-path
difference, not streaming; making `runPca` consume the DelayedArray directly is
an open change and is what the 1.3M rung actually needs.

Mixed read paths make a cost comparison meaningless. Pick one for the ladder.

## 7. GPU accounting

* `gpu_upload` / `gpu_download` are charged to **compute**, not I/O. PCIe
  transfer is part of using the GPU; bucketing it with disk reads would let a
  GPU arm look cheaper on compute than it is.
* **VRAM is not captured at all.** `performance.txt` reports host RSS only, so
  for a rapids arm the memory figure is an undercount of unknown size. The
  number you want is VRAM + host. This needs `nvidia-smi` sampling alongside
  the job; until it exists, no memory comparison involving a GPU arm is
  like-for-like.

## 8. Do not assume a single power law

Over two or more orders of magnitude you will cross memory cliffs and index
ceilings (the R-side 2^31 limit bites at the top TENx rungs). Fit segmented and
report breakpoints; a line drawn through a discontinuity produces a slope that
describes nothing, and R² will not always catch it.

## 9. denet is not wired up

The `-prof` entrypoints (`pca-prof: denet pca.py`, etc.) exist in every module
and **have never run**. Three things block them:

1. ob invokes `./{entrypoint}`, so `denet pca.py` becomes `./denet pca.py` —
   it looks for a file called `denet` inside the module directory.
2. denet is not in any conda env; it is at `~/.cargo/bin/denet`.
3. denet 0.6.0 requires a subcommand: `denet run <cmd>`, not `denet <cmd>`.
   `denet /bin/echo hi` returns "unrecognized subcommand".

So all phase timings so far come from `obkit-events.jsonl` (start/end events the
modules emit) plus `performance.txt`, not from denet. denet would add sampled
thread activity and memory pressure at ~100ms, which is exactly what section 5
had to infer — worth fixing, but it is three fixes.

## 10. The stage-attribution trap (hit three times)

The output tree nests `STAGE/module/.hash` all the way down, so a CLUST node's
path **contains** `/PCA/` and `/NNG/`. Any `re.search(r"/(PCA|NNG|CLUST)/")` or
`"/PCA/" in path` matches the wrong stage and fails silently:

* `fit_scaling.py` attributed every NNG and CLUST event to PCA;
* `pairwise.py` read the clustering seed as the PCA seed, which would have made
  the bootstrap resample the nested unit and report bands far too narrow;
* `collect.py` needed the same fix for ancestor parameters.

**Always take the last match, or match on structure** (`anc.parent.parent.name
== "PCA"`). And `DATA_dataset_name` is null for omni-data modules, which
silently dropped 696 tm-droplet nodes from every group_by — use the `dataset`
column `collect.py` now emits.

## 11. A protocol that would be publishable

1. 6–8 rungs spanning ≥2 orders of magnitude (5k → 500k+).
2. 3 repeats per rung.
3. Serial, idle machine, cold cache, fixed governor.
4. One read path for every arm.
5. Fit time and RAM separately, per phase, with CIs; segmented if the residuals
   show structure.
6. Verify effective cores per arm from `performance.txt` before interpreting.
7. State n\* for any predicted cost, and keep it inside the measured span.
