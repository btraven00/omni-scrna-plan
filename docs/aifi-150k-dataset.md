# A 150k-cell arm from the AIFI immune atlas

Why this document exists: every cost number in the clustering slice today is
measured at 12k–55k cells, and at that scale the cost axis is measuring startup,
not algorithms. Fitting job CPU time against cell count across the three current
datasets gives log-log slopes of

    nn-scanpy    0.01      cl-scrapper  0.17      pc-scanpy    0.38
    cl-rapids    0.17      cl-scanpy    0.38      cl-seurat    0.44
    pc-rapids   -0.59

where 1.0 would be linear in cells. `nn-scanpy` is the most expensive stage in
the chain at ~17s and it is **completely flat** from 12k to 55k: 4.5x the data
costs it nothing. `pc-rapids` is *negative*, which is only possible if the
measurement is dominated by CUDA context setup that a bigger input amortises
better; its peak RSS sits at 1.47-1.57 GB regardless of n, i.e. the allocator
pool, not the data.

So the current conclusion "rapids costs more and saves no time" is really
"rapids pays a fixed GPU init that 55k cells never amortises". Whether that
survives at 150k+ is an open question, and it is the main reason for this arm.

Two compounding limits on the same numbers, both worth remembering when reading
any cost figure from the existing runs:

* **Everything ran single-threaded.** `slices/clustering.yaml` declares no
  `resources: cores:` on any stage, and measured `cpu_time / wall` is below 1.0
  for every module (cl-seurat 0.92, nn-scanpy 0.57, fi-scrapper 0.43). The
  `--cores 4` on the command line is snakemake job concurrency, not threads.
* **Counter noise is large relative to the gaps.** Repeat runs of the *same*
  CLUST config differ by p90/p10 = 1.17x in time (CV 6.6%), while the largest
  time ratio between any two pipelines is 1.39x. Memory is the opposite: RSS
  repeatability is under 1%.

## The source

`/home/b/phd/data/4015a993-cac9-4732-94ff-74c48b37dc29.h5ad` — an AIFI human
immune atlas, CELLxGENE export.

    1,234,234 cells x 32,357 genes
    X = CSR float32, all-integer  -> RAW UMI COUNTS
    2.04e9 nonzeros (~1,652 per cell)
    no layers, no .raw, obsm/X_umap present

**Use this file, not `4015a993.hvg2k.h5ad`.** The hvg2k derivative is already
cut to 2000 HVGs, so FEAT would be selecting 2000 out of 2000 -- a no-op that
silently removes a stage from the comparison while still emitting a
`normalized_selected_h5` that looks correct.

### Labels

Three nested levels plus an ontology column, and no rare-class problem at L1/L2:

    AIFI_L1       9 classes   smallest  537 cells   none under 100
    AIFI_L2      29 classes   smallest  345 cells   none under 100
    AIFI_L3      71 classes   smallest   53 cells   one under 100
    cell_type    38 classes   CELLxGENE ontology terms

The hierarchy is the real asset. Running L2 and L3 on *identical cells* varies
label granularity with everything else held fixed, which separates "does the
method scale" from "does the label set get harder". No other dataset here can do
that -- tm-facs has one label set at 81 classes and that is all.

### Natural subsets

`donor_id` has 10 levels; `cohort`, `assay` and `disease` are all single-level,
so donor is the only structure available.

    BR1025 165,105   BR1041 141,078   BR1035  34,702
    BR1030 159,980   BR1048 140,917   BR1033  34,170
    BR1014 155,477   BR1008  99,095
    BR1037 152,654   BR1010 151,056

Take **BR1010, 151,056 cells**, for the first arm: it is the size the cost
question needs, it removes the batch confound, and it leaves the other donors as
future rungs (two ~ 320k, the full atlas 1.23M).

**Do NOT select "the smallest donor"** -- the distribution is far wider than the
top of it suggests, and the smallest is 34,170 cells, below tm-facs and useless
for the scale question. Name the donor explicitly.

## Build plan

Adapt `build_pbmc_omnidata.py`, which already writes the exact four-file DATA
contract omni-data serves from a `file://` directory:

    {id}.h5ad                  counts
    {id}.clusters_truth.tsv    cell_id \t truths
    {id}.clusters_truth_num.txt
    {id}_properties.yaml       batch_var / sample_var / labels_var

Two differences from the pbmc build: the source is an h5ad (no `scx` step), and
the subsetting must **stream**. 2.04e9 nonzeros will not load; read `indptr`
once (10 MB), then pull each selected row's slice and rebuild CSR incrementally.
The output is ~248M nonzeros, about 2.0 GB.

Three things that fail silently if skipped:

1. **Write `obs` flat.** anndata 0.13 emits `nullable-string-array` and
   categorical encodings that anndataR *drops with only a WARNING*, so the R
   modules (`fi-scrapper`, `cl-seurat`, `cl-metrics-r`) would see missing columns
   and the run would complete while producing nonsense. Plain fixed-length
   strings.
2. **Drop blank/NA labels before counting**, or `clusters_truth_num.txt` is
   wrong and every downstream metric is scored against the wrong truth count.
3. **Decide about `var/feature_is_filtered`.** CELLxGENE zeroes those genes in
   X. Leaving the zero columns in changes what "2000 HVGs" is selected *from*.
   Drop at build time or handle in FILT, but do it deliberately.

## Wiring into the slice

Same shape as the tm-facs DATA entry:

    - id: aifi-d1-l2
      software_environment: omni_data
      repository:
        url: https://github.com/btraven00/omni-data
        commit: 80611e83d13c8bd0d15fbc033c6a5e1a7900a0b8
      parameters:
        - uri: file:///home/b/phd/data/aifi_omnidata/aifi-d1-l2/
          uri_type: directory

`file://` is fine for development and **not** for results: it is not
reproducible, and changing the uri later changes the DATA parameter hash and
discards the entire computed tree below it. Publish to
omnibenchmark.mls.uzh.ch before anything from this arm is reported.

## What will hurt

* **Drop `GRAPH-M` (`nn-metrics-r`) for this arm.** It is in the slice at line
  378 and already runs 63s at 55k; graph metrics scale badly and this is the
  most likely stage to make the run unaffordable.
* **`CLUST-M` is the throughput risk**, not any single job -- it is ~76% of the
  job count in a sweep of this shape.
* **`FILT` peaked at 5.9 GB on the 44k arm.** Budget well above that.
* **Declare `resources: cores:` per stage.** There is real headroom unused (see
  above), though Leiden itself will not benefit.

## Trim the sweep first

The current grid is 20 PCA seeds x 5 clustering seeds x 3 resolutions x 4 arms.
That is not a first run at 150k. Start with **3 PCA x 3 clustering seeds at the
three complete resolutions** to get cost and a sanity check, then decide whether
the full seed grid is affordable.

Keep the **crossing**. A reduced crossed design still decomposes into variance
components; an OFAT one does not, and the whole
[variance decomposition](../../omni-scrna-analysis/plot_variance.py) depends on
it. Note also that only resolutions carrying the *complete* seed grid are usable
downstream -- the analysis scripts filter on that, because a 3-seed exploratory
cell wins a "best median resolution" contest on noise.

## What this arm does NOT fix

These labels come from annotation that itself used clustering, so this arm is
**circular in exactly the way tm-facs was added to avoid**. It extends the size
axis; it does not replace tm-facs, and the two answer different questions. It is
also one tissue and one assay, so it adds no tissue diversity either.

Related: [`scalability-fitting.md`](scalability-fitting.md) for why slopes are
the characteristic quantity and intercepts are machine-dependent, and
[`nondominated-pipelines.md`](nondominated-pipelines.md) for how the cost axis
enters the two-stage selection.

See also [`duo2018-calibration-tier.md`](duo2018-calibration-tier.md): this arm
is the *scale* tier of three, and the tier table there says which question each
dataset is allowed to answer.
