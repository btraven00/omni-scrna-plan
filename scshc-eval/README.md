# scSHC cluster-significance eval

Self-contained probe of [scSHC](https://github.com/igrabski/sc-SHC)'s `testClusters`
before wiring it into the metrics module. It takes raw counts + a clustering and
reports how many clusters are statistically significant (non-significant splits
get merged).

## Contents
- `data/data.h5ad` — real dataset counts (29606 cells × 36753 genes, raw in `layers/counts`)
- `data/…cl-scanpy.tsv` — a real clustering from the out_ci run (scanpy, K=13, 26830 cells)
- `scshc_eval.R` — worker: loads counts once, runs `testClusters` per TSV
- `run.sh` — runs the bundled TSV(s); pass extra TSV paths to sweep many
- `environment.yml` — conda/micromamba env (`r-scshc` from the almost-conductor channel)

## Run
```sh
rsync -a scshc-eval/ server:scshc-eval/          # ~820 MB (the h5ad)
ssh server
cd scshc-eval
conda env create -f environment.yml && conda activate scshc-eval   # once
SCSHC_CORES=16 ./run.sh                                             # full dataset
```
The full dataset (~27k cells, K=13) is slow — expect many minutes per TSV, hence
the server. `run.sh` defaults to every `data/*.tsv`; pass explicit paths to pick.

## Output (`out/`)
- `summary.tsv` — one row per TSV: `id, n_cells, input_K, surviving_K, retained_frac, elapsed_s, status` (appended live, survives crashes)
- `<id>.json` — per-TSV detail: the `merges` (which input clusters collapsed) and the printed significance `tree` (node names carry the split p-values)

Read-off: `surviving_K < input_K` ⇒ some clusters weren't statistically distinct.
On the fixture all 8 truth cell lines survive (K 8→8) — it's clean, so nothing
merges. To exercise the merge path, feed an over-clustered TSV.

## Running on real omnibenchmark output
The cluster TSVs are `data_clusters.tsv` under the `CLUST/` nesting, and the counts
are the `data.h5ad` at the dataset root. Point the worker straight at them:
```sh
ROOT=/path/to/out_ci/DATA/data/.<hash>
SCSHC_CORES=16 Rscript scshc_eval.R "$ROOT/data.h5ad" out \
  $(find "$ROOT" -name data_clusters.tsv)
```
Heads-up: the full dataset (~27k cells) is far slower per TSV than the fixture —
that's why this runs on a server.
