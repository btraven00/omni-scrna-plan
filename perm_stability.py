#!/usr/bin/env python3
"""Pairwise-ARI stability across cell-order permutations, per clustering arm.

The stability variable is ARI BETWEEN permutations, not ARI to truth: it needs
no label set, so it separates "is this method order-stable" from "is this
method accurate" -- two things a single ARI-to-truth number conflates. Accuracy
is reported alongside, never instead.

Works straight off the cluster TSVs, so it does not wait on CLUST-M: cl-sc3s
sits on CLUST-E and the resolver still binds cl-metrics-r to CLUST only, which
would otherwise leave the sc3s arm with no metrics at all.

Every TSV is keyed by cell_id, so arms are joined on barcode rather than on row
position -- which is the whole point, since row position is what varies.

    python perm_stability.py out_pbmc
    python perm_stability.py out_pbmc pbmc.h5ad celltype.l1 celltype.l2 protein_cluster

A dataset can carry more than one truth (pbmc has three annotation levels plus a
protein-derived one), and stability and accuracy answer different questions --
so every truth given is scored, next to the one stability number.
"""
import csv
import os
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from sklearn.metrics import adjusted_rand_score as ari

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "out_pbmc")


def load(p):
    # stdlib csv, not polars: labels stay strings, which ARI handles fine
    with open(p) as f:
        r = csv.reader(f, delimiter="\t")
        next(r)
        return {row[0]: row[-1] for row in r}


def h5ad_obs(path, col):
    """One obs column of an .h5ad as {barcode: label}, categorical or not."""
    import h5py
    dec = lambda b: b.decode() if isinstance(b, bytes) else str(b)
    with h5py.File(path) as f:
        ids = f["obs"][f["obs"].attrs["_index"]][:]
        g = f["obs"][col]
        if isinstance(g, h5py.Group):          # categorical: codes -> categories
            cats = g["categories"][:]
            vals = [cats[c] for c in g["codes"][:]]
        else:
            vals = g[:]
    return {dec(i): dec(v) for i, v in zip(ids, vals)}


def truths():
    """{name: {barcode: label}} -- the DATA stage's own truth, plus any asked for."""
    out = {}
    t = next(OUT.rglob("*.clusters_truth.tsv"), None)
    if t:
        out["(DATA truth)"] = load(t)
    for col in sys.argv[3:]:                   # argv[2] = the .h5ad they live in
        out[col] = h5ad_obs(sys.argv[2], col)
    return out


TRUTHS = truths()


def arms(out):
    """Map arm -> {permutation_seed: labels-by-barcode}.

    The graph arms carry permutation_seed on the NNG directory; sc3s carries it
    on its own CLUST-E directory. Two ob axes, one experiment -- joined here by
    seed VALUE, which is exactly the join the collector would have to do.
    """
    # ob stores results in hashed dirs (.02c668ab) and exposes the parameters
    # only as SIBLING SYMLINKS (..._permutation_seed-7_random_seed-42). So the
    # seed is never in a real path -- walk in through the symlinks or find
    # nothing.
    found = {}
    for link in out.rglob("*permutation_seed-*"):
        if not link.is_symlink():
            continue
        seed = int(re.search(r"permutation_seed-(\d+)", link.name).group(1))
        # Everything in the symlink name EXCEPT permutation_seed is the config.
        # It has to key the arm, or the four n_runs values (and the two Leiden
        # resolutions) collapse into one ECDF -- pooling a curve index with the
        # random variable, which is exactly the mistake the design forbids.
        cfg = re.sub(r"_?permutation_seed-\d+", "", link.name)
        real = link.parent / os.readlink(link)
        for tsv in real.rglob("*_clusters.tsv"):
            s = str(tsv)
            stage = "CLUST-E" if "/CLUST-E/" in s else "CLUST"
            mod = re.search(rf"/{stage}/([^/]+)/", s)
            if mod:
                found.setdefault(f"{mod.group(1)}  [{cfg}]", {})[seed] = tsv
    return found


for arm, by_seed in sorted(arms(OUT).items()):
    seeds = sorted(by_seed)
    if len(seeds) < 2:
        continue
    lab = {s: load(by_seed[s]) for s in seeds}
    bc = sorted(set.intersection(*(set(v) for v in lab.values())))
    vec = {s: np.array([lab[s][b] for b in bc]) for s in seeds}
    ks = {s: len(set(vec[s])) for s in seeds}

    pw = [ari(vec[a], vec[b]) for a, b in combinations(seeds, 2)]
    print(f"\n=== {arm}  ({len(seeds)} permutations, {len(bc)} cells)")
    print(f"  pairwise ARI   mean {np.mean(pw):.4f}  sd {np.std(pw):.4f}  "
          f"min {np.min(pw):.4f}  max {np.max(pw):.4f}   ({len(pw)} pairs)")
    if 0 in vec:   # seed 0 is the identity control, not a privileged answer
        v0 = [ari(vec[0], vec[s]) for s in seeds if s]
        print(f"  vs identity    mean {np.mean(v0):.4f}  min {np.min(v0):.4f}")
    kv = sorted(set(ks.values()))
    print(f"  k              {kv}  moved {sum(k != ks[seeds[0]] for k in ks.values())}"
          f"/{len(seeds)}")

    # Accuracy, reported next to stability and never instead of it: a method
    # that always returns the same wrong partition is perfectly stable. Computed
    # here rather than read from CLUST-M because cl-sc3s sits on CLUST-E and the
    # resolver does not feed it to cl-metrics-r.
    for name, t in TRUTHS.items():
        # a truth file with no barcode overlap (the ladder has no labels) would
        # otherwise print a triumphant ARI 1.0000 over zero cells
        keep = [b for b in bc if b in t]
        if not keep:
            continue
        tv = np.array([t[b] for b in keep])
        acc = [ari(tv, np.array([lab[s][b] for b in keep])) for s in seeds]
        print(f"  vs {name:<16s} ARI mean {np.mean(acc):.4f}  sd {np.std(acc):.4f}  "
              f"spread {np.max(acc) - np.min(acc):.4f}  "
              f"(k_truth {len(set(tv))}, {len(keep)} cells)")

    # Cost comes free: ob writes performance.txt beside every module output.
    perf = [p for s in seeds
            for p in [by_seed[s].parent / f"{by_seed[s].name.split('_clusters')[0]}_performance.txt"]
            if p.exists()]
    if perf:
        rows = [dict(zip(*(l.split("\t") for l in p.read_text().splitlines()[:2])))
                for p in perf]
        w = [float(r["s"]) for r in rows]; m = [float(r["max_rss"]) for r in rows]
        print(f"  cost           {np.mean(w):.1f}s/run  {np.mean(m):.0f}MB max_rss")
