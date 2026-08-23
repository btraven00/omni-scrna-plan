#!/usr/bin/env python3
"""kNN seed sweep in one process: read the PCA once, loop seeds in memory.

The benchmark slice (benchmark_knnseed.yaml) answers the same question but pays
per seed for a conda activation, a 149MB TSV parse, a 53MB HDF5 write and an R
process. This does the loop in-process and keeps nothing on disk but the result
table -- ~4x faster and no 6GB of intermediate graphs.

Numbers are only comparable within one environment: the same seed run under
pynndescent 0.5.13 and 0.6.0 disagrees on ~18% of kNN edges, which is more than
the seed itself moves (~5%). Do not compare a sweep run here against the
benchmark's -- run it in the benchmark's own conda env (which has leidenalg,
so --leiden igraph) if you need to line the two up.

Also records the things that plausibly change the answer across machines
(platform, ISA, numba threading layer, BLAS) into the output header, so the
same script run on linux-x86_64 and macOS-arm64 can be diffed directly.

  ./knn_seed_sweep.py --pcas_tsv X_pcas.tsv --truth X.clusters_truth.tsv \
      --seeds 1-100 --out sweep.tsv
"""

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pcas_tsv", type=Path, help="PCA embedding TSV (cell_id + PCs)")
    p.add_argument("--truth", type=Path, help="truth TSV: cell_id <TAB> label")
    p.add_argument("--out", type=Path, default=Path("knn_seed_sweep.tsv"))
    p.add_argument("--seeds", default="1-100", help="'1-100' or '1,7,42'")
    p.add_argument("--n_neighbors", type=int, default=15)
    p.add_argument("--flavor", default="umap", choices=["umap", "gauss"])
    p.add_argument("--transformer", default="auto",
                   choices=["auto", "pynndescent", "sklearn"],
                   help="'auto' = scanpy's own choice: exact under 8192 cells, "
                        "approximate above -- so 'auto' on a small input gives "
                        "100 identical rows, by design")
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--leiden", default="igraph", choices=["igraph", "rapids"],
                   help="'rapids' needs rapids_singlecell and a GPU")
    p.add_argument("--leiden_seed", type=int, default=42,
                   help="held fixed so the kNN seed is the only thing moving")
    p.add_argument("--n_jobs", type=int, default=1,
                   help="scanpy settings.n_jobs; pynndescent stays reproducible "
                        "at a fixed seed for n_jobs > 1")
    p.add_argument("--selftest", action="store_true")
    return p.parse_args()


def expand_seeds(spec):
    if "-" in spec and "," not in spec:
        lo, hi = spec.split("-")
        return list(range(int(lo), int(hi) + 1))
    return [int(s) for s in spec.split(",")]


def environment():
    """What would make two machines disagree. Recorded, not controlled."""
    import numba
    import scanpy as sc
    import sklearn
    env = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "sklearn": sklearn.__version__,
        "numba": numba.__version__,
    }
    try:
        import pynndescent
        env["pynndescent"] = pynndescent.__version__
    except ImportError:
        env["pynndescent"] = None
    try:
        from numba.np.ufunc import parallel
        env["numba_threading_layer"] = parallel.threading_layer()
    except Exception:
        # only resolves after a parallel region has actually run
        env["numba_threading_layer"] = "unresolved"
    try:
        env["blas"] = [d["internal_api"] + "-" + d["version"]
                       for d in np.show_config(mode="dicts")["Build Dependencies"]
                       .get("blas", {}).get("name", "")] or None
    except Exception:
        env["blas"] = np.show_config(mode="dicts").get(
            "Build Dependencies", {}).get("blas", {}).get("name")
    return env


def cluster(adata, how, resolution, seed):
    if how == "rapids":
        import rapids_singlecell as rsc
        rsc.get.anndata_to_GPU(adata)
        rsc.tl.leiden(adata, resolution=resolution, obsp="connectivities",
                      random_state=seed, key_added="cluster")
        rsc.get.anndata_to_CPU(adata, convert_all=True)
    else:
        import scanpy as sc
        sc.tl.leiden(adata, resolution=resolution, random_state=seed,
                     flavor="igraph", n_iterations=2, directed=False,
                     key_added="cluster")
    return adata.obs["cluster"].astype(str).to_numpy()


def sweep(embedding, cell_ids, truth, args):
    import anndata as ad
    import scanpy as sc
    from sklearn.metrics import (adjusted_mutual_info_score, adjusted_rand_score,
                                 v_measure_score)

    sc.settings.n_jobs = args.n_jobs
    sc.settings.verbosity = 0
    transformer = None if args.transformer == "auto" else args.transformer

    rows = []
    for i, seed in enumerate(args.seeds, 1):
        t0 = time.perf_counter()
        adata = ad.AnnData(X=np.zeros((embedding.shape[0], 1)))
        adata.obs_names = cell_ids
        adata.obsm["X_pca"] = embedding
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, method=args.flavor,
                        use_rep="X_pca", random_state=seed, transformer=transformer)
        labels = cluster(adata, args.leiden, args.resolution, args.leiden_seed)
        rows.append({
            "seed": seed,
            "k_found": len(np.unique(labels)),
            "ARI": adjusted_rand_score(truth, labels),
            "AMI": adjusted_mutual_info_score(truth, labels),
            "VM": v_measure_score(truth, labels),
            "secs": round(time.perf_counter() - t0, 2),
        })
        print(f"  [{i}/{len(args.seeds)}] seed {seed}: ARI {rows[-1]['ARI']:.4f} "
              f"k={rows[-1]['k_found']} ({rows[-1]['secs']}s)", flush=True)
    return pd.DataFrame(rows)


def summarize(df):
    from scipy import stats
    a = df["ARI"].to_numpy()
    if len(a) < 3:
        return {"n": len(a)}
    return {
        "n": len(a), "mean": a.mean(), "sd": a.std(ddof=1),
        "min": a.min(), "max": a.max(), "spread": a.max() - a.min(),
        "cv_pct": a.std(ddof=1) / a.mean() * 100,
        "skew": float(stats.skew(a)), "excess_kurtosis": float(stats.kurtosis(a)),
        "shapiro_p": float(stats.shapiro(a).pvalue),
    }


def selftest():
    """Same seed twice must give the same ARI; different seeds may not."""
    import types
    rng = np.random.default_rng(0)
    n = 900
    emb = np.vstack([rng.normal(c, 1.0, (n // 3, 8)) for c in (-4, 0, 4)])
    truth = np.repeat(["a", "b", "c"], n // 3)
    ids = [f"c{i}" for i in range(n)]
    args = types.SimpleNamespace(
        seeds=[7, 7, 8], n_neighbors=15, flavor="umap", transformer="pynndescent",
        resolution=1.0, leiden="igraph", leiden_seed=42, n_jobs=1)
    df = sweep(emb, ids, truth, args)
    assert df.ARI[0] == df.ARI[1], f"seed 7 not reproducible: {df.ARI[0]} vs {df.ARI[1]}"
    assert df.ARI.max() > 0.9, f"3 separated blobs should be easy, got {df.ARI.max()}"
    print("selftest ok")


def main():
    args = parse_args()
    if args.selftest:
        return selftest()
    if not (args.pcas_tsv and args.truth):
        sys.exit("--pcas_tsv and --truth are required (or use --selftest)")
    args.seeds = expand_seeds(args.seeds)

    pcas = pd.read_csv(args.pcas_tsv, sep="\t", index_col=0)
    truth_df = pd.read_csv(args.truth, sep="\t", index_col=0)
    truth = truth_df.iloc[:, 0].reindex(pcas.index)
    if truth.isna().any():
        sys.exit(f"{int(truth.isna().sum())} cells in the PCA have no truth label")
    print(f"{pcas.shape[0]} cells x {pcas.shape[1]} PCs, "
          f"{truth.nunique()} truth labels, {len(args.seeds)} seeds")

    env = environment()
    print("env: " + json.dumps(env))
    df = sweep(pcas.to_numpy(dtype=np.float64), list(pcas.index),
               truth.to_numpy(), args)

    summary = summarize(df)
    with open(args.out, "w") as fh:
        fh.write("# env: " + json.dumps(env) + "\n")
        fh.write("# args: " + json.dumps({k: v for k, v in vars(args).items()
                                          if k != "seeds"}, default=str) + "\n")
        fh.write("# summary: " + json.dumps(summary, default=float) + "\n")
        df.to_csv(fh, sep="\t", index=False)
    print(f"wrote {args.out}")
    print("  " + "  ".join(f"{k}={v:.4f}" if isinstance(v, float) else f"{k}={v}"
                           for k, v in summary.items()))


if __name__ == "__main__":
    main()
