# Duò 2018 as a calibration tier

[Duò, Robinson & Soneson 2018](https://doi.org/10.12688/f1000research.15666.3),
*A systematic performance evaluation of clustering methods for single-cell
RNA-seq data*. Same group as this benchmark's listed benchmarker, so an arm on
these datasets is directly comparable to published reference numbers in a way
tm-facs and the AIFI atlas are not.

Distributed as the `DuoClustering2018` Bioconductor ExperimentHub package.

## The datasets, measured

Enumerated locally 2026-09-20 from `DuoClustering2018` 1.18.0, not recalled from
the paper:

    dataset                  cells   features  classes  min class
    Koh                        531     48,981        9         22
    KohTCC                     531    811,938        9         22
    Kumar                      246     45,159        3         69
    KumarTCC                   263    803,405        3         83
    SimKumar4easy              500     43,606        4         38
    SimKumar4hard              499     43,638        4         61
    SimKumar8hard              499     43,601        8         16
    Trapnell                   222     41,111        3         72
    TrapnellTCC                227    684,953        3         73
    Zhengmix4eq              3,994     15,568        4        997
    Zhengmix4uneq            6,498     16,443        4        500
    Zhengmix8eq              3,994     15,716        8        398

## The blocker: almost all of this is below the exact-kNN line

scanpy builds an **exact** neighbour graph below 8,192 cells (and the ANN seed is
inert below that), so on these datasets the graph is a constant and the kNN seed
axis does not exist. Worse for our purposes, the 4,096-cell figure recorded in
[the divergence budget](../../omni-knn-seed) is the point below which this has
been verified byte-identical across seeds.

**Only `Zhengmix4uneq` (6,498) clears 4,096.** `Zhengmix4eq` and `Zhengmix8eq`
sit at 3,994 — *102 cells* below the line. Everything else is in the hundreds.

So the variance decomposition in
[`plot_variance.py`](../../omni-scrna-analysis/plot_variance.py) would be
**degenerate by construction** on this tier, not merely small: one of its two
factors cannot move. Likewise every cost conclusion — these are already smaller
than pancreas (12,115), which is itself in the regime where measured scaling
slopes are ~0 and the cost axis is reading startup.

**This tier is for accuracy and calibration only.** Exclude it explicitly from
cost and seed-sensitivity analysis rather than letting it quietly dilute them.

## Two further traps

* **The four `*TCC` variants are not gene matrices.** They are transcript
  compatibility equivalence classes: 685k–812k features. NORM and FEAT would be
  operating on a completely different object. Exclude them unless the TCC
  representation is itself the question.
* **The `filteredExpr10` / `filteredHVG10` / `filteredM3Drop10` variants have
  gene selection already applied** — the same defect as `4015a993.hvg2k.h5ad`
  (see [`aifi-150k-dataset.md`](aifi-150k-dataset.md)). FEAT would become a
  no-op that still emits a plausible-looking `normalized_selected_h5`. Use the
  `sce_full_*` objects.

## What is genuinely worth having

**`SimKumar4easy` vs `SimKumar4hard` is a clean difficulty contrast.** 500 vs
499 cells, 4 classes both, ~43.6k genes both. The *only* thing that differs is
group separation. That is a difficulty axis with nothing else confounded, and
nothing else in this benchmark provides one.

Note `SimKumar8hard` does **not** extend that cleanly: it doubles the classes at
fixed n, so its minimum class falls to 16 cells. It confounds difficulty with a
rare-class problem. Report it separately from the 4easy/4hard pair or not at all.

**The truth labels are not circular.** SimKumar\* labels are simulation ground
truth. Zhengmix\* are FACS-purified populations mixed computationally, so the
labels come from protein sorting, not from clustering the same expression the
methods will cluster. Both are *stronger* on this axis than tm-facs, whose
labels are expert annotation. For the chance→oracle bracket figure, which asks
"does the embedding carry the labels at all", non-circular truth is exactly what
the argument needs.

Caveat on the small classes: Koh has a 22-cell class and SimKumar8hard a
16-cell one. Partition metrics are unstable on classes that small, and some
consensus methods drop clusters under 20 cells outright.

## The three tiers

Each answers a different question; none substitutes for another.

| tier | datasets | n | truth | used for |
|---|---|---|---|---|
| calibration | Duò `sce_full_*` (non-TCC) | 222–6,498 | simulation / FACS | accuracy, difficulty gradient, comparability to published numbers |
| real, non-circular | tm-facs | 44,104 | expert annotation | accuracy and stability on real data |
| scale | AIFI donor subset | ~151,056 | annotation (circular) | cost, scaling slopes, where rapids may overcome |

pancreas and tm-droplet remain as the incumbent accuracy arms but share the
circularity of the AIFI tier.

## Environment note

`bioconductor-*` conda packages ship the R tarball under `share/` and rely on a
post-link script; `AnnotationHub` additionally needs **`r-dbplyr < 2.4`** (2.4+
raises `Arguments in ... must be used / ..1 = Inf`). Installing
`bioconductor-singlecellexperiment` incrementally into a solved env broke the R
binary (`libgfortran.so.3`) — solve everything in one go. Working spec:

    micromamba create -p ENV --channel-priority flexible -c conda-forge -c bioconda \
      "bioconductor-duoclustering2018=1.18.0=r43hdfd78af_0" \
      "bioconductor-singlecellexperiment" "r-dbplyr=2.3.4"

Pin this into `envs/` before the tier becomes a real arm.
