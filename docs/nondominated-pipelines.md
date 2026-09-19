# Non-dominated pipelines: accuracy, stability, cost

Extends [`order-stability-dominance.md`](order-stability-dominance.md), which
defines the ECDF, the cluster bootstrap and the ε-violation statistic for a
*single* axis. This one is about combining axes: which pipelines survive once
accuracy, stability, scalability and resources are all on the table.

Read that document first. Everything it says about ε, about resampling runs
rather than pairs, and about stability not being accuracy, holds here unchanged.

## 1. The unit

A **pipeline** is a PCA → NNG → CLUST path with its parameters. Pre-processing
(DATA/FILT/NORM/FEAT) is deliberately *not* part of it:

* it is shared by every pipeline in a sweep, so per-pipeline attribution is
  arbitrary — in the run measured below, 2 PCA nodes feed 272 CLUST nodes;
* its cost is dominated by I/O churn rather than by the method under test;
* including it would swamp the axis that discriminates. Measured here: FEAT
  peaks at **16.8 GB** against PCA's 0.36 GB, so a pipeline "memory" number
  including FEAT would be the same number for every pipeline.

Cost is therefore summed over PCA + NNG + CLUST only.

## 2. The axes

| axis | estimand | direction | units |
|---|---|---|---|
| accuracy | median ARI over the randomisation draws | higher | fraction of the chance→oracle range |
| stability | pairwise ARI *between* runs (§3 below) | higher | ARI |
| scalability | slope of log(time) on log(n_cells) | lower | dimensionless |
| time | median wall seconds, PCA+NNG+CLUST | lower | s |
| memory | max RSS over those three stages | lower | GB |

**Accuracy is normalised by the CLUSTBOUND bounds**, not reported raw. The
chance floor and the supervised ceiling turn ARI into "fraction of the
achievable range claimed", which is what makes a tolerance expressible: a 0.003
ARI difference is 0.3% of range on pancreas and meaningless, and without the
denominator nothing says so. The bounds also give the ECDF its left and right
anchors.

## 3. Stability is between runs, not dispersion of accuracy

The estimand is the one in the sibling document: ARI(π_i, π_j) for two runs of
the same configuration under independently drawn randomisation. It needs no
ground truth, so it separates *reproducible* from *correct*.

> **A pipeline that is reproducibly wrong must not score well.** Measured here:
> pancreas at resolution 1.0 has modal-k stability 0.99–1.00 across 20 seeds and
> ARI 0.27 against an oracle of 0.97. On a naive (accuracy, stability) front it
> is non-dominated *on stability* and would be recommended. Gate stability
> behind an accuracy floor, or read it only in oracle-normalised units.

`pareto.py` currently reports **accuracy dispersion** (p90−p10 of ARI over
seeds) as a stand-in. That is not the same quantity and cannot make the
separation above; it is labelled as a proxy in the script and should be replaced
by `cluster/pairwise_stability.py`'s pairwise ARI once it is wired to the parquet.

## 4. Sources of variation to pool

The ECDF pools only *genuine randomisation*: sources that are exchangeable draws
from a population that exists.

* **pool**: clustering seed, PCA seed, cell subsample (80%, without replacement)
* **do not pool**: resolution, k, n_neighbors — these are design choices, not
  draws. Pooling them makes the ECDF a property of the grid's density rather
  than of the algorithm; doubling the points below 0.05 would move it with no
  change to the method. Facet by them instead.

Keep the design **crossed and balanced**. 20 subsamples × 3 PCA seeds and
3 × 20 pool to different distributions at the same total n, so a ragged design
re-introduces the weighting problem by the back door — and a crossed design is
what lets the variance be decomposed per knob afterwards.

Resolution belongs in the report as an **envelope**, not a distribution:
best-over-grid (what tuning buys) and worst-over-grid (what you get without it).
Measured on be1-subset, that envelope is the one method-level difference in this
data that clears the effect floor:

| arm | worst over grid | best over grid | width |
|---|---|---|---|
| cl-rapids | 0.5937 | 0.9702 | 0.377 |
| cl-scanpy | 0.6147 | 0.9660 | 0.351 |
| cl-scrapper | 0.6421 | 0.9573 | 0.315 |

scrapper has the lowest ceiling and the highest floor: the flattest response to
resolution, i.e. the least tuning-sensitive.

## 5. Tolerances

With 4–5 axes almost everything is non-dominated, so ε-dominance is not a
refinement — it is what makes the front a decision instead of a list.

| axis | tolerance | form | rationale |
|---|---|---|---|
| accuracy | 5% of chance→oracle range | absolute | a perfectly consistent 0.3%-of-range win is not a result |
| stability | 0.05 ARI | absolute | matches the ε≲0.05 band of the sibling doc |
| time | 1.5× | multiplicative | below this, scheduling noise and contention dominate |
| memory | 10× | multiplicative | only catastrophic differences change a decision |

State the asymmetry rather than hiding it in the numbers: **10× on memory means
memory almost never binds.** That is defensible on a workstation and wrong on a
shared cluster with per-job limits; the multiplier is the place to encode which
you mean.

**Measured: the memory multiplier is the only tolerance that changes the answer
at this scale.** Time never binds (max observed ratio 1.36×, inside either 1.5×
or 2×). Memory spans 3.03×, so:

| tolerance | non-dominated |
|---|---|
| mem ×10, time ×2 | 8 — rapids on the front for both datasets |
| mem ×2, time ×2 | 6 — rapids drops out of both |

`cl-rapids` peaks at 1.48 GB against 0.49 GB for scanpy/scrapper and buys
+0.008 accuracy, which is *inside* the 0.05 accuracy tolerance and therefore a
tie. Its place on the ×10 front was an artifact of the multiplier, not a
property of the method. Use **×2 at this scale**.

> **`max_rss` is host memory only.** The GPU arm's VRAM is not captured by
> `performance.txt`, so `cl-rapids` is simultaneously understated (no VRAM) and
> already 3× the others on host RSS. Any memory comparison involving a GPU arm
> needs `nvidia-smi` sampling or it compares different quantities.

## 6. What this data can and cannot measure

Run: `slices/clustering.yaml`, 2 datasets × 3 clusterers × 2 resolutions × 20
clustering seeds.

**Measurable now** — accuracy, accuracy-dispersion, time, memory. The front:

```
Non-dominated (acc ±0.05, disp ±0.05, time ×1.5, mem ×10)
  be1-subset  res 0.05  cl-rapids    acc 0.964  disp 0.005  17.8s  1.48 GB
  be1-subset  res 0.05  cl-scrapper  acc 0.961  disp 0.005  16.9s  0.49 GB
  be1-subset  res 0.05  cl-scanpy    acc 0.956  disp 0.005  15.3s  0.49 GB
  pancreas    res 0.05  cl-rapids    acc 0.578  disp 0.047  35.8s  1.47 GB
  pancreas    res 0.05  cl-scanpy    acc 0.577  disp 0.006  33.2s  0.57 GB
  pancreas    res 0.05  cl-scrapper  acc 0.577  disp 0.006  35.0s  0.57 GB
```

Every resolution-1.0 configuration is dominated, and within resolution 0.05 all
three clusterers survive: they differ by <1% of range in accuracy, <1.2× in
time and <3× in memory, all inside tolerance. **The honest reading is that the
resolution choice is the decision and the clusterer choice is not** — accuracy
spans 0.27→0.58 across resolutions on pancreas and 0.001 across clusterers.

**Not measurable here, and why**

* **Scalability slope** — needs ≥3 cell counts per pipeline. This run has one
  size per dataset. The TENx ladder is the instrument; note that some arms
  truncate (R's 2³¹ index ceiling at the top rungs, rapids needs a GPU), and a
  slope from 3 rungs is not comparable to one from 6. Report slope *and*
  intercept: a good slope with a bad constant loses at every size actually run.
* **Subsample and PCA-seed variation** — neither is swept yet. The FILT module
  has no `--subsample_frac`; adding it is the one piece of module code this
  design needs, and everything downstream inherits the draw. PCA seed is a
  parameter already but frozen at one value here.
* **Stability proper** — §3; needs the partitions, not the parquet.

## 7. Caveats to write down now

* **The timings in §6 were collected on a contended machine** (`--cores 4`, four
  jobs in flight). The sibling document's requirement — the machine otherwise
  idle — was not met. They are good enough to show that time does not bind at
  this scale; they are not publication numbers. Re-measure serially before any
  cost claim.
* **Two different ε are in play.** `order-stability-dominance.md` uses
  Leshno–Levy on **CDF** differences, ∫_V (F_A − F_B) dt ⁄ ∫ |F_A − F_B| dt.
  `dominance.py` currently implements Dror et al. on squared **quantile**
  differences. Both are "almost first-order dominance" and they do not agree
  numerically. The CDF form is the house convention; the script should move to
  it, and until it does the two must not be compared or reported under one name.
* **ε is magnitude-weighted, not rank-based** (in either form). Nineteen narrow
  losses against one large win reads as dominance. Bounding the score does not
  fix it: one ARI of 1.0 against nineteen at 0.50 beats a constant 0.55 at
  ε ≈ 0.43, i.e. inside the conventional τ = 0.5. Prefer a stricter τ and read
  the ECDF beside the scalar. Pinned in `test_dominance.py`.
* **All arms here share one kNN graph** (`nn-scanpy`), so method differences are
  implementation, not graph. A fair clusterer comparison needs each on its own
  graph, or the shared graph stated as the design.
* **n = 20 seeds** resolves a proportion to about ±0.05. Enough to rank, too
  coarse to separate 0.95 from 0.98.

## 12. Adding the self-ensembling arms (sc3s, Markov stability)

The point of the framework is this comparison: methods that buy stability by
spending compute, priced against sweeping Leiden harder. Both modules exist.

### Where they attach

| arm | stage | consumes | note |
|---|---|---|---|
| sc3s | CLUST-E | `embedding_tsv` | reads the embedding, not the graph, so it is CLUST-E not CLUST. `permutation_seed` goes here, not on NNG. |
| pygenstability (`cl-pgs`) | CLUST | `neighbors_h5` | Markov stability over scales; the scale selector is argmin-NVI |

Both land in CLUST-M through the shared `clusters_tsv` fan-in, so no new metric
stage is needed — the same poem call scores them beside Leiden and beside the
CLUSTBOUND floor/ceiling.

### Fix k and n_components; do not let the method choose

sc3s can pick k itself. **Turn that off** and give it the same k the comparison
is held at. Two reasons:

* k is sc3s' INPUT and Leiden's OUTPUT. A method whose k never moves scores
  trivially well on any k-stability measure, and that is not a merit — it is the
  one axis on which the arms are not comparable at all.
* Letting it choose introduces a second free parameter that differs per arm,
  and the resolution<->k mapping is already scale-dependent per dataset.

Same for **`n_components`: hold it at the value every other arm uses (50 here,
or 30 — pick one and state it).** How PCA truncation affects clustering is a
real question and a separate experiment; varying it inside this comparison
would confound the ensembling question with the truncation question and answer
neither.

### Cost accounting

Use `--cost ensemble` for the headline and `--cost sweep` beside it (see §11 of
elite.py's docstring). sc3s at `n_runs=N` and pygenstability over scales pay
their ensemble internally in one job; Leiden pays its sweep across k jobs and
its seed spread across n more. Reporting only one mode picks the winner by
accounting convention.

The figure that answers the question: **stability gain against compute
increase** — delta pairwise-ARI on the y axis, delta (seconds, peak RSS) on the
x, with single-run Leiden at the origin. An ensembling method earns its place
only if it sits above the line Leiden traces as its own ensemble grows.

### Expect a negative result, and design for it

Prior measurements in this project, worth re-testing rather than assuming:

* consensus (sc3s) measured **less** order-stable than Leiden at matched k;
* the Markov-stability arm was **dead as a stability arm** — its same-graph
  control (0.70) was as bad as its permutation spread (0.66), so the argmin-NVI
  scale selector was what moved, not Louvain.

"You pay 10x compute for no stability gain" is a publishable finding. Set the
run up to measure it, not to confirm the opposite.

### Practical blockers

* scSHC is blocked by **n^2 memory**, not time — it will not reach the larger
  datasets.
* `testClusters` silently returns K=1 when the root split isolates a tiny
  cluster; drop clusters below ~20 cells before scoring.
