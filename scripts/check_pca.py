#!/usr/bin/env python3
"""Check PCA h5 outputs for validity and compare subspace similarity across runs."""

import subprocess
import re
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import csv

import h5py
import numpy as np
from scipy.spatial import procrustes

EXPECTED_DATASETS = ["cell_ids", "embedding", "gene_ids", "loadings", "variance", "variance_ratio"]


def load_embedding(path: str, standardize: bool = False) -> np.ndarray:
    """Return embedding as (cells, components), normalizing orientation.

    standardize=True column-standardizes (zero mean, unit variance per PC)
    before returning — required for cross-module Procrustes because
    OnlinePCA.jl uses a per-PC scale convention (λ/M) that differs from
    scanpy/scrapper (σ), so each component has its own scale factor that a
    single global Procrustes scalar cannot absorb.
    """
    with h5py.File(path, "r") as f:
        emb = f["embedding"][:]
        n_cells = f["cell_ids"].shape[0]
    if emb.shape[0] != n_cells:
        emb = emb.T
    if standardize:
        emb = (emb - emb.mean(axis=0)) / (emb.std(axis=0) + 1e-12)
    return emb


def check_file(path: str) -> list[str]:
    issues = []
    try:
        with h5py.File(path, "r") as f:
            missing = [k for k in EXPECTED_DATASETS if k not in f]
            if missing:
                issues.append(f"missing keys: {missing}")
                return issues

            embedding = f["embedding"][:]
            loadings = f["loadings"][:]
            variance = f["variance"][:]
            variance_ratio = f["variance_ratio"][:]
            n_components = variance.shape[0]
            n_cells = f["cell_ids"].shape[0]
            n_genes = f["gene_ids"].shape[0]

            if set(embedding.shape) != {n_components, n_cells}:
                issues.append(f"embedding shape {embedding.shape} inconsistent with n_components={n_components}, n_cells={n_cells}")
            if set(loadings.shape) != {n_components, n_genes}:
                issues.append(f"loadings shape {loadings.shape} inconsistent with n_components={n_components}, n_genes={n_genes}")

            for name, arr in [("embedding", embedding), ("loadings", loadings),
                               ("variance", variance), ("variance_ratio", variance_ratio)]:
                if np.any(np.isnan(arr)):
                    issues.append(f"{name} contains NaN")
                if np.any(np.isinf(arr)):
                    issues.append(f"{name} contains Inf")

            if not np.all(variance >= 0):
                issues.append("variance has negative values")
            if not np.all(variance_ratio >= 0):
                issues.append("variance_ratio has negative values")
            if not np.all(np.diff(variance) <= 0):
                issues.append("variance is not monotonically decreasing")

    except Exception as e:
        issues.append(f"could not open file: {e}")

    return issues


def label(path: str) -> str:
    name = Path(path).name
    m = re.match(r"pca-(\w+)_(.+)\.h5", name)
    if m:
        return f"{m.group(1)}/{m.group(2)}"
    return name


def upstream_key(path: str) -> str:
    m = re.match(r"(.+/five-pca[^/]*)/", path)
    return m.group(1) if m else path


def load_timing(h5_path: str) -> dict | None:
    perf = Path(h5_path).parent / "datasets_performance.txt"
    if not perf.exists():
        return None
    with open(perf) as f:
        reader = csv.DictReader(f, delimiter="\t")
        row = next(reader, None)
    if row is None:
        return None
    return {"wall_s": float(row["s"]), "cpu_s": float(row["cpu_time"])}


def main():
    result = subprocess.run(
        ["find", "out", "-name", "pca-*.h5"],
        capture_output=True, text=True,
    )
    files = sorted(result.stdout.strip().splitlines())

    if not files:
        print("No pca-*.h5 files found under out/")
        raise SystemExit(1)

    # --- validity checks ---
    print("=== validity ===")
    ok = failed = 0
    for path in files:
        issues = check_file(path)
        status = "OK" if not issues else "FAIL"
        print(f"[{status}] {path}")
        for issue in issues:
            print(f"       ! {issue}")
        if issues:
            failed += 1
        else:
            ok += 1
    print(f"\n{ok} OK, {failed} failed out of {len(files)} files")

    # --- timings ---
    print("\n=== timings (wall seconds) ===")
    print(f"  {'method/solver':<38} {'mean':>7} {'std':>7}  n")
    print("─" * 60)
    by_solver: dict[str, list[float]] = defaultdict(list)
    for path in files:
        t = load_timing(path)
        if t is None:
            continue
        by_solver[label(path)].append(t["wall_s"])
    for lbl, vals in sorted(by_solver.items()):
        print(f"  {lbl:<38} {np.mean(vals):>6.1f}s {np.std(vals):>6.1f}s  {len(vals)}")

    # --- procrustes comparison ---
    print("\n=== procrustes disparity (0=identical, 1=orthogonal) ===")
    groups: dict[str, list[str]] = defaultdict(list)
    for f in files:
        groups[upstream_key(f)].append(f)

    # collect all disparities per pair type across upstream groups
    pair_disparities: dict[str, list[float]] = defaultdict(list)
    for key, paths in sorted(groups.items()):
        if len(paths) < 2:
            continue
        for a, b in combinations(paths, 2):
            try:
                emb_a = load_embedding(a, standardize=True)
                emb_b = load_embedding(b, standardize=True)
                if emb_a.shape != emb_b.shape:
                    continue
                _, _, disparity = procrustes(emb_a, emb_b)
                pair_key = " vs ".join(sorted([label(a), label(b)]))
                pair_disparities[pair_key].append(disparity)
            except Exception:
                pass

    print(f"  {'pair':<55} {'mean':>7} {'std':>7}  n")
    print("─" * 78)
    for pair, vals in sorted(pair_disparities.items()):
        print(f"  {pair:<55} {np.mean(vals):>7.4f} {np.std(vals):>7.4f}  {len(vals)}")

    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
