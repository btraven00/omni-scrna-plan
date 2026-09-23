# Comparing Order-Stability by Stochastic Dominance

How to compare clustering methods on their sensitivity to cell storage order,
when each method has a different knob and a different compute cost.

Written for the seeds slice (`benchmark_seeds.yaml`), which sweeps
`permutation_seed` over 20 values with every other seed frozen. The arms are
`cl-scanpy` / `cl-rapids` (Leiden, `permutation_seed` on NNG) and `cl-sc3s`
(consensus k-means, `permutation_seed` on CLUST-E because it reads `pcas_tsv`
rather than the graph).

## The estimand

For arm *A* in configuration *c*:

> **S**<sub>A,c</sub> = ARI(π<sub>i</sub>, π<sub>j</sub>), where π<sub>i</sub> and π<sub>j</sub> are the partitions *A*
> produces under two independently drawn cell orderings.

The target is its CDF, **F**<sub>A,c</sub>(t) = P(S ≤ t).

Two properties make this the right variable rather than ARI-against-truth:

* It needs **no ground truth**, so it separates *is this method order-stable*
  from *is this method accurate*. A single ARI-to-truth number conflates them,
  and a method that reliably returns the same wrong partition scores 1.0 here
  and badly on accuracy. Report accuracy alongside, never instead.
* It **privileges no ordering**. Comparing the 19 permuted runs against
  `permutation_seed-0` treats the unpermuted order as the reference answer; it
  is not one, it is one arbitrary draw. (On full pbmc, seed 0 lands mid
  distribution.) All 190 pairs use every ordering symmetrically.

Higher ARI is better, so **lower CDF = more stable**, and *A dominates B* means
*A's curve lies entirely below B's*. Put that convention on the figure — it is
the easiest thing for a reader to invert.

## 1. The ECDF

F̂(t) = (1/190) · #{ i<j : ARI<sub>ij</sub> ≤ t }, over the 20·19/2 = 190 pairs.

Consistent for F. The point estimate is fine; only the *uncertainty* breaks if
the 190 values are treated as independent (see §2).

Plot:

| element | choice | why |
|---|---|---|
| curves | 7 — three Leiden arms, four sc3s `n_runs` | the sc3s family marching right IS the frontier |
| x range | the occupied range (~0.45–0.95), not [0,1] | the crossing region is narrow and is the whole point |
| colour | by family; lightness by `n_runs` | shows compute buying stability at a glance |
| annotation | cost (s/run) on every curve | cost indexes the comparison, so it belongs on the figure |
| caption | state n = 20 | 190 steps look smooth and invite over-reading |

## 2. Cluster bootstrap

**The trap.** The 190 pairwise ARIs come from only 20 underlying runs — each run
appears in 19 pairs. Effective sample size is ~20, not 190. Resampling *pairs*
shatters the dependence and produces bands far too narrow.

**The fix.** Resample the independent units, which are the permutations:

```python
M = pairwise_ari_matrix(partitions)        # 20x20, computed ONCE
for b in range(B):
    idx = rng.choice(20, size=20, replace=True)
    pairs = [M[a, c] for a, c in combinations(idx, 2) if a != c]
    F_b = ecdf(pairs)
```

Three details that are not cosmetic:

* **Precompute `M`.** Recomputing ARI inside the loop is 190 × B evaluations on
  156,881 cells — hours. Indexing a cached 20×20 matrix makes the whole
  bootstrap microseconds, so B = 2000 costs nothing.
* **Drop self-pairs.** Sampling with replacement puts the same permutation in
  `idx` twice, and ARI(π, π) = 1.0 exactly. Keeping those biases every replicate
  toward stability. That is what `if a != c` is for.
* **Pointwise bands are not simultaneous.** Percentile bands at each *t* are 190
  pointwise intervals, not a joint region, and cannot carry a dominance claim.
  Draw them labelled *pointwise*, and put the inference on the scalar ε below,
  where the bootstrap distribution is one-dimensional and a percentile CI is
  honest.

## 3. The epsilon-violation statistic (almost-FSD)

Strict first-order dominance requires F<sub>A</sub>(t) ≤ F<sub>B</sub>(t) for **all** *t*.

> **Overlapping supports do NOT imply a crossing.** It is tempting to conclude
> from sc3s' max (0.7101) exceeding cl-scanpy's min (0.6948) that FSD must fail.
> It does not follow, and here it is false: at t = 0.70, F<sub>scanpy</sub> ≈ 0.01 while
> F<sub>sc3s</sub> ≈ 0.98 — the curves are nowhere near each other, and measured ε is
> exactly 0. Overlap is necessary for a crossing, not sufficient. Compute ε;
> do not predict it from the supports.

The binary FSD test is still the wrong instrument when curves *do* approach,
because it discards how small the violation is. Leshno–Levy's ε quantifies it:

> **V** = { t : F<sub>A</sub>(t) > F<sub>B</sub>(t) }
>
> **ε** = ∫<sub>V</sub> (F<sub>A</sub> − F<sub>B</sub>) dt ⁄ ∫ |F<sub>A</sub> − F<sub>B</sub>| dt

The area where dominance fails, as a fraction of the total area between the
curves.

| ε | reading |
|---|---|
| 0 | strict FSD holds |
| ≲ 0.05 | almost first-order dominance; violation confined to a negligible region |
| → 0.5 | the curves genuinely cross; no dominance statement available |

Both curves are step functions, so integration is **exact and trivial**: build
the common grid from the sorted union of all observed ARI values, evaluate both
curves, sum the rectangles between consecutive knots. No quadrature.

Compute ε in **both directions** and report the smaller with its direction —
ε(A over B) and ε(B over A) are not complements.

For inference, compute ε inside the same bootstrap loop of §2. The claim becomes

> cl-scanpy almost-dominates sc3s at `n_runs=5`: ε = 0.02, 95% CI [0.00, 0.06]

which is a number with uncertainty, rather than a failed binary test.

## 4. What to report

Four comparisons — each sc3s configuration against the best Leiden arm
(`cl-scanpy`). Not all 21 pairs: that is multiplicity for no gain.

| arm | config | cost s/run | mean pairwise ARI | ε vs cl-scanpy | 95% CI |
|---|---|---|---|---|---|

Plus the figure, plus the headline the sweep exists to produce:

> the `n_runs` at which ε stops being small, expressed as a multiple of a single
> Leiden pass.

Cost is **not** a random variable and does not get an ECDF. It is deterministic
per configuration, so it *indexes* the family of dominance verdicts. That family
is the frontier. `performance.txt` beside every module output already carries
wall seconds and `max_rss`, so the cost axis needs no instrumentation — but the
machine must be otherwise idle when the timings are collected.

## 5. Caveats to write down now, not discover later

* **Leshno–Levy ε has a known critique** (Tzeng, Huang & Shih 2013): the
  threshold does not map cleanly onto a well-behaved utility class. Use ε
  descriptively — "the violation region is 2% of the area between the curves" —
  and claim no decision-theoretic force. That reading is unimpeachable and
  sufficient.
* **Stability dominance is not accuracy dominance.** Run the same machinery on
  ARI-vs-truth separately. sc3s may escape stability dominance at high `n_runs`
  while remaining accuracy-dominated.
* **The k confound is real and matched-k reverses the reading.** Instability
  rises with k in both families, so scoring sc3s' k=20 against Leiden at k≈30
  compares sc3s against Leiden at its worst. Measured on full pbmc: both Leiden
  arms get *more* stable at k=20 than at their default resolution (cl-scanpy
  0.7780 → 0.8299, cl-rapids 0.7249 → 0.8282). Correcting for it widened the
  gap from 0.04 to 0.09 and removed the only contested comparison in the table.
  A matched-k arm costs 8s per job. Always include one.

* **The resolution ↔ k mapping is SCALE-DEPENDENT. Never carry it between
  datasets.** Calibrate on the target dataset's own graph, at one seed, before
  adding a matched-k arm. Measured here:

  | resolution | k at 20,000 cells | k at 156,881 cells |
  |---|---|---|
  | 0.4 | 13 | **20** |
  | 0.8 | 16 | 27 |
  | 1.0 | 19 | 29 |
  | 1.2 | **20** | ~35 |

  Resolution 1.2 is the matched-k arm at 20k and is nowhere near it at 156k —
  reading the wrong column puts the "matched" arm *further* from the target k
  than the default arm it replaced, while looking like a correction. This
  mistake was made here and caught only by calibrating. One `cluster` call per
  candidate resolution on an existing `neighbors.h5` is enough; validate the
  environment first by checking it reproduces a k you already know.
* **sc3s' `k` never moves, and that is not a merit.** k is its input, not its
  output. That row is trivially perfect and is the one axis on which the arms
  are not comparable at all. Only Leiden and scSHC choose their own k.
* **n = 20 gives a 20-step ECDF.** Enough for gross dominance, thin exactly in
  the tails where the crossing happens. If the crossing point carries the
  result, add permutations to the two arms that bracket it rather than to all.

## 6. Self-checks to build in

* ε(A, A) = 0 by construction.
* ε = 0 for synthetically separated supports.
* Bootstrap CI = [0, 0] for an arm against itself.
