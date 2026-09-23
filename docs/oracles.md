# The CLUSTBOUND Bounds: a Chance Floor and a Supervised Ceiling

What the two bound modules compute, why the ceiling is swept over two
classifiers, and how to read the result.

An ARI of 0.63 means nothing on its own. It could be a weak clusterer on a
clean embedding, or the best any method could do on an embedding that does not
carry the labels. CLUSTBOUND answers that by bracketing every arm from below
and above, on the *same* embedding the arm was given.

## They emit partitions, not scores

Both modules write `{dataset}_clusters.tsv` — `cell_id<TAB>cluster` — under the
same `clusters_tsv` output id that CLUST emits. They score nothing themselves.
CLUST-M then picks them up and runs the identical poem call it runs on a real
clustering arm.

This is the load-bearing design choice. A bound and an arm are the same number
*by construction* rather than because two implementations agree — and they did
not agree: sklearn normalises AMI by the arithmetic mean of the entropies,
poem by the max. A ceiling computed with sklearn and arms scored with poem
would have been quietly on different scales.

The cost is a wiring requirement: CLUST and CLUSTBOUND both produce
`clusters_tsv`, so CLUST-M needs an `ob` with shared-output-id fan-in
(`gather-lean-3`, after commit `b78abe9`). On `main`, CLUST-M silently takes
one producer and the bounds vanish with no error.

## The floor — `cb-chance`

A uniform random partition at *k*:

```python
k = n_clusters if n_clusters > 0 else len(np.unique(y))
rng.integers(0, k, size=len(cell_ids))
```

ARI and AMI are chance-corrected, so this scores ≈0 by construction. The point
is not the centre, it is the **spread** — the band a weak arm has to clear.
One run is one draw, so the distribution comes from the seed sweep
(`random_seed: [1..5]`), the same way `permutation_seed` is swept on NNG.

`n_clusters: [0, 10, 20]` is swept because chance level depends on *k* — the
same confounder that makes ARI at a fixed resolution rank arms by their cluster
count.

**Measured, on the definitive run:** the floor is not binding at these sizes.

| dataset | median ARI | max ARI | sd |
|---|---|---|---|
| aifi-d1 | −0.00000 | 0.00001 | 0.00001 |
| pancreas | 0.00000 | 0.00034 | 0.00014 |
| tm-droplet | −0.00001 | 0.00001 | 0.00002 |
| tm-facs | −0.00002 | 0.00006 | 0.00004 |

With n ≥ 12k and k ≤ 55, chance ARI is zero to four decimal places; pancreas
shows the widest spread because it is the smallest dataset. The floor earns its
place as a *check* — if an arm ever landed in that band, something is broken —
not as a threshold anyone has to clear.

**What it deliberately is not:** a gene-permutation null (shuffle each gene
across cells, keep the marginals) is the stronger test, because it asks whether
the pipeline manufactures structure out of none. But it has to re-enter the DAG
at FILT and rerun everything downstream, which makes it a benchmark arm, not a
module on this stage.

## The ceiling — `cb-oracle`

Cross-validated predictions from a classifier that was handed the truth labels:

```python
cv   = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
pred = cross_val_predict(clf, X, y, cv=cv)
```

`X` is the embedding, `y` the truth labels. The predicted labels are written out
as if they were a clustering, and scored as one.

It is **not a strict upper bound**. An unsupervised method could in principle
beat a given classifier. What it answers is narrower and more useful: *is the
clusterer weak, or does the embedding simply not carry the labels?* If the
oracle is also low, no arm was ever going to score well and the fault is
upstream of CLUST — in FEAT, NORM, or the data.

Classes with fewer than `n_splits` members are dropped before the CV rather
than letting StratifiedKFold raise; single-digit classes are normal here (tm-facs
carries 81 labels, several tiny). CLUST-M inner-joins on `cell_id` and reports
the shortfall as `n_dropped`.

## Why two classifiers

Because the ceiling is only ever as good as the classifier that defines it, and
a classifier that fails on a given embedding produces a ceiling that is not just
imprecise but **actively misleading** — it reads as "the embedding does not
carry the labels" about an embedding that plainly does. That is the one error
this stage must not make.

The measurement that settled it, on tm-facs-marrow (5,015 cells, 50 PCs, 22
labels):

| classifier | ARI |
|---|---|
| logreg | 0.88 |
| knn | 0.67 |
| HistGradientBoosting | **0.07** |

Boosting does not move between 100 and 400 rounds — it converges somewhere
degenerate rather than underfitting. A 0.07 ceiling sits *below* every real arm.
So `hgb` remains available but is not swept, and the two that are swept are
chosen to fail differently: **logreg** is linear and global, **knn** (k=15) is
local and non-parametric. An embedding where cell types are linearly separable
favours the first; one where they are curved but locally compact favours the
second.

The rule is therefore **take the max over classifiers, never one arm**.

**What the definitive run actually shows** (median ARI over 84 embeddings per
cell):

| dataset | knn | logreg | margin |
|---|---|---|---|
| aifi-d1 | 0.7319 | **0.8524** | 0.121 |
| pancreas | 0.9583 | **0.9723** | 0.014 |
| tm-droplet | 0.9549 | **0.9677** | 0.013 |
| tm-facs | 0.8897 | **0.8917** | 0.002 |

logreg wins on all four, so in hindsight the second arm bought insurance rather
than a better number. The insurance came closest to mattering on aifi-d1, where
using knn alone would have understated the ceiling by 0.12 — and on tm-facs,
where the two agree to 0.002 but knn recovers only 79 of the 81 labels against
logreg's 81, losing two classes outright.

Two arms is also cheap: 84 embeddings × 2 classifiers × 5 folds, against the
5,290 metric rows each dataset carries in total.

## The ceiling is a property of the data, not of the PCA

Bounds are computed per embedding, so each one is itself a distribution over the
PCA grid — 20 seeds × 3 PCA modules per dataset. That distribution is very
tight:

| dataset | min | median | max | range |
|---|---|---|---|---|
| aifi-d1 | 0.8513 | 0.8524 | 0.8547 | 0.0034 |
| pancreas | 0.9703 | 0.9723 | 0.9726 | 0.0023 |
| tm-droplet | 0.9673 | 0.9677 | 0.9685 | 0.0012 |
| tm-facs | 0.8912 | 0.8917 | 0.8931 | 0.0019 |

Across 84 embeddings the ceiling moves by less than 0.004. Which PCA you used
does not meaningfully change what is recoverable — so a single ceiling per
dataset is a fair summary, and any large movement in an arm's ARI across
embeddings is the arm's, not the embedding's.

## Reading the bracket

The headroom, not the raw ARI, is the number worth reporting:

| dataset | best arm | ceiling | % of ceiling |
|---|---|---|---|
| aifi-d1 | 0.736 | 0.852 | 86% |
| pancreas | 0.633 | 0.972 | **65%** |
| tm-droplet | 0.844 | 0.968 | 87% |
| tm-facs | 0.742 | 0.892 | 83% |

Pancreas is the case that justifies the whole stage. Its best unsupervised arm
scores 0.63, the lowest of the four — which on its own reads as the hardest
dataset. The ceiling says the opposite: the embedding carries the labels almost
perfectly (0.97, the highest of the four), and a third of what is there is being
left on the table by the clusterers. Without the ceiling that is invisible; with
it, pancreas is where method improvement has the most to gain.

## Configuration

```yaml
- id: cb-chance
  parameters:
    - bound: chance
      n_clusters: [0, 10, 20]        # 0 = the truth's own k
      random_seed: [1, 2, 3, 4, 5]   # the sweep IS the distribution

- id: cb-oracle
  parameters:
    - bound: oracle
      classifier: [logreg, knn]      # take the max; see above
      n_splits: 5
      random_seed: 42
```

Both run from `bounds/bounds.py` in `split-stages-plan-modules/metrics`
(entrypoint `bounds-py`), stage inputs `embedding_tsv` and
`rawdata_clusters_truth`.

## Caveat

**AMI is unusable on aifi-d1** — 3,705 of its rows are NaN, all four other
datasets are clean. Use ARI for anything involving aifi, including these
brackets, until that is tracked down.
