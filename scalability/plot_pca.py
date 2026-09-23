#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas", "matplotlib"]
# ///
"""Plot the PCA scalability curve from an `ob collect performance` table.

    ob collect performance out_run                 # writes out_run/performances.tsv
    ./scalability/plot_pca.py out_run/performances.tsv

x is cells (parsed from the tenx-NNNNk rung id), y is wall time and peak RSS.
One series per module x solver; upstream stages are dropped -- PCA is the only
thing under test here.
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

tsv = sys.argv[1] if len(sys.argv) > 1 else "out_run/performances.tsv"
out = sys.argv[2] if len(sys.argv) > 2 else "pca_scalability.png"

df = pd.read_csv(tsv, sep="\t")
df = df[df["stage"] == "PCA"].copy()
df["cells"] = df["dataset"].str.extract(r"(\d+)k", expand=False).astype(int) * 1000
df["solver"] = df["params"].map(lambda p: json.loads(p)["PCA"]["solver"])
df["series"] = df["module"] + ":" + df["solver"]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
for name, g in df.groupby("series"):
    g = g.sort_values("cells")
    axes[0].plot(g["cells"], g["s"], "o-", label=name)
    axes[1].plot(g["cells"], g["max_rss"], "o-", label=name)
for ax, ylab in zip(axes, ["wall time (s)", "peak RSS (MB)"]):
    ax.set(xscale="log", yscale="log", xlabel="cells", ylabel=ylab)
    ax.grid(alpha=.3, which="both")
axes[0].legend(fontsize=8)
fig.suptitle("PCA scalability: TENx ladder")
fig.savefig(out, dpi=150)
print(f"wrote {out}  ({len(df)} PCA runs, {df['cells'].nunique()} rungs)")
