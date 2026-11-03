#!/usr/bin/env python
"""Fail unless the .h5 holds a sparse TENx matrix as scrapper writes it.

Reference is HDF5Array::writeTENxMatrix (scrapper, seurat): matrix/{data float64,
indices uint32, indptr int64, shape uint32, genes, barcodes}, genes x cells CSC.
Dtype widths may differ (scanpy writes int32/vlen strings) as long as indptr can
index every stored value. Storage only: dense values written sparse pass.

Scale note (5M-cell runs): nr-scanpy (normalize.py) casts indptr with
.astype(np.uint32), which wraps silently past 4.29e9 stored values; the indptr
check below is what catches it. fe-scanpy (feat-select.py) writes scipy's
indptr as-is, which scipy widens to int64 by itself. The fix for both is
indptr.astype(np.int64), matching scrapper.
"""
import sys

import h5py
import numpy as np

KEYS = {"data", "indices", "indptr", "shape", "genes", "barcodes"}
# what writeTENxMatrix writes; anything else is a WARN, not a FAIL
REF = {"data": "float64", "indices": "uint32", "indptr": "int64", "shape": "uint32"}


def validate_file(path):
    try:
        with h5py.File(path, "r") as f:
            g = f["matrix"]
            if not isinstance(g, h5py.Group):
                raise ValueError("matrix is a dense dataset, expected a sparse group")
            missing = KEYS - set(g)
            if missing:
                raise ValueError(f"matrix/ is missing {', '.join(sorted(missing))}")
            if g["data"].dtype.kind != "f":
                raise ValueError(f"data is {g['data'].dtype}, expected float")
            for k in ("indices", "indptr", "shape"):
                if g[k].dtype.kind not in "iu":
                    raise ValueError(f"{k} is {g[k].dtype}, expected integer")
            for k in ("genes", "barcodes"):
                if g[k].dtype.kind not in "SO":
                    raise ValueError(f"{k} is {g[k].dtype}, expected strings")
            n_genes, n_cells = g["shape"][:]
            nnz = g["data"].shape[0]
            if g["indptr"].shape[0] != n_cells + 1:
                raise ValueError(f"indptr has {g['indptr'].shape[0]} entries, expected {n_cells + 1}")
            if g["indices"].shape[0] != nnz:
                raise ValueError("data and indices differ in length")
            if np.iinfo(g["indptr"].dtype).max < nnz:
                raise ValueError(f"indptr {g['indptr'].dtype} cannot index {nnz} values")
            if g["genes"].shape[0] != n_genes or g["barcodes"].shape[0] != n_cells:
                raise ValueError("genes/barcodes length does not match shape")
            for k, want in REF.items():
                if g[k].dtype != want:
                    print(f"WARN\t{path}\t{k} is {g[k].dtype}, scrapper writes {want}")
            for k in ("genes", "barcodes"):
                if g[k].dtype.kind != "S":
                    print(f"WARN\t{path}\t{k} is variable-length, scrapper writes fixed bytes")
    except Exception as e:
        print(f"FAIL\t{path}\t{e}")
        return
    print(f"OK\t{path}")


for arg in sys.argv[1:]:
    validate_file(arg)
