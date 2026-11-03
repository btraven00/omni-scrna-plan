#!/usr/bin/env python
"""Validate a one-data stage .h5ad.

The downstream R modules read obs columns (batch/sample/labels) with
rhdf5::h5read, which expects a flat array. anndata writes string columns as
categorical *groups* (categories + codes), which h5read mangles into a
degenerate vector. Fail if any obs column is categorical-encoded.
"""
import sys

import h5py


def validate_file(path):
    try:
        with h5py.File(path, "r") as f:
            bad = [c for c in f["obs"]
                   if f["obs"][c].attrs.get("encoding-type") == "categorical"]
            if bad:
                raise ValueError(
                    "obs columns are categorical-encoded (must be flat "
                    f"string-arrays for rhdf5::h5read): {', '.join(bad)}")
    except Exception as e:
        print(f"FAIL\t{path}\t{e}")
        return
    print(f"OK\t{path}")


for arg in sys.argv[1:]:
    validate_file(arg)
