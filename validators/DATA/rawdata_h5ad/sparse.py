#!/usr/bin/env python
"""Fail if a DATA stage .h5ad stores layers["counts"] dense."""
import sys

import scanpy as sc
from scipy.sparse import issparse

for path in sys.argv[1:]:
    # TODO: loads the whole file; read lazily if larger inputs make this slow
    layers = sc.read_h5ad(path).layers
    if "counts" not in layers:
        print(f"FAIL\t{path}\tno layers['counts']")
    elif not issparse(layers["counts"]):
        print(f"FAIL\t{path}\tlayers['counts'] is {type(layers['counts']).__name__}, expected sparse")
    else:
        print(f"OK\t{path}")
