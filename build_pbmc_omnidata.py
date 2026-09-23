#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["h5py"]
# ///
"""Recreate the pbmc omni-data DATA dir from the Azimuth h5seurat.

pbmc_multimodal.h5seurat --(scx, SCT/counts)--> pbmc_omnidata/
    pbmc.h5ad                 DATA contract: counts in layers/counts (CSR), no X
    pbmc.clusters_truth.tsv   cell_id \t <celltype.l2>
    pbmc.clusters_truth_num.txt
    pbmc_properties.yaml

NB: SCT/counts are SCTransform-*corrected* counts (depth-regressed), not raw UMIs.
Fine to drive the pipeline; for real significance read-offs use GEO GSE164378 raw.

Run (uv auto-installs h5py into an ephemeral env):
    uv run build_pbmc_omnidata.py        # or ./build_pbmc_omnidata.py
Still needs `scx` and `curl` on PATH.
"""
import csv
import hashlib
import os
import subprocess

import h5py

URL = "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat"
# Pin after the first download: set to the sha256 fetch_src() prints. Left None
# = trust-on-first-use (the atlas ships no checksum of its own).
SRC_SHA256 = "ff154ce3672a33b91b5b05297147554aae69887ba670cd42af6e27a94682bbaa"
SRC = "/home/b/phd/data/pbmc_multimodal.h5seurat"
SCX = "/home/b/phd/scx/target/release/scx"
OUT = "/home/b/phd/data/pbmc_omnidata"
LABEL = "celltype.l2"  # -> clusters_truth and labels_var
PROPS = {"batch_var": "donor", "sample_var": "orig.ident", "labels_var": LABEL}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def add_obs_column(f, name, values):
    """Add/replace a flat string-array obs column and register it in column-order
    (so scanpy/anndata sees it as a column, not stray data)."""
    sd = h5py.string_dtype(encoding="utf-8")
    obs = f["obs"]
    if name in obs:
        del obs[name]
    ds = obs.create_dataset(name, data=list(values), dtype=sd)
    ds.attrs["encoding-type"] = "string-array"
    ds.attrs["encoding-version"] = "0.2.0"
    co = obs.attrs.get("column-order")
    order = [c.decode() if isinstance(c, bytes) else str(c)
             for c in (list(co) if co is not None else [])]
    if name not in order:
        order.append(name)
    if "column-order" in obs.attrs:
        del obs.attrs["column-order"]
    obs.attrs.create("column-order", order, dtype=sd)


def fetch_src():
    """Download the Azimuth PBMC reference if absent. Resumes into a .part file
    so an interrupted download never masquerades as a complete SRC, then verifies
    the sha256 against SRC_SHA256 (or prints it to pin on first run)."""
    if not os.path.exists(SRC):
        os.makedirs(os.path.dirname(SRC), exist_ok=True)
        part = SRC + ".part"
        print(f"fetching {URL}")
        subprocess.run(["curl", "-fL", "-C", "-", "-o", part, URL], check=True)
        os.replace(part, SRC)
    got = sha256(SRC)
    if SRC_SHA256 is None:
        print(f"source sha256 (pin this in SRC_SHA256): {got}")
    elif got != SRC_SHA256:
        raise SystemExit(f"source checksum mismatch:\n  got      {got}\n  expected {SRC_SHA256}")


def main():
    fetch_src()
    os.makedirs(OUT, exist_ok=True)
    h5ad = os.path.join(OUT, "pbmc.h5ad")

    # 1. scx: SCT corrected counts -> h5ad. --exclude drops every slot the DATA
    #    stage doesn't consume at conversion time (skips reading obsm too, and
    #    never writes the two 161k^2 obsp graphs), leaving just X + obs + var.
    #    f32 (default) stores these small integer counts exactly; the validator's
    #    floor check passes. Needs the scx slot-filter build (feature/convert-
    #    slot-filter or later).
    subprocess.run(
        [SCX, "convert", SRC, h5ad, "--assay", "SCT", "--layer", "counts",
         "--exclude", "obsm,varm,obsp,varp,uns,layers"],
        check=True,
    )

    # 2. The one thing --exclude can't do: scx writes counts to X, the DATA
    #    contract wants them in layers/counts with no top-level X. Move it.
    #    (scx already emits obs/_index + var/_index and string-array obs cols.)
    with h5py.File(h5ad, "r+") as f:
        layers = f.require_group("layers")
        layers.attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
        f.move("X", "layers/counts")
        assert "X" not in f and "layers/counts" in f, "counts not routed to layer"
        assert "_index" in f["obs"] and "_index" in f["var"], "index not named _index"
        cell_ids = f["obs/_index"].asstr()[:]
        labels = f[f"obs/{LABEL}"].asstr()[:]

        # Optional: merge protein-modality clusters (run cluster_protein.R first)
        # as an ADT-based ground-truth obs column, aligned to obs order by barcode.
        prot_tsv = os.path.join(OUT, "pbmc.protein_clusters.tsv")
        if os.path.exists(prot_tsv):
            with open(prot_tsv) as fh:
                m = {r["cell_id"]: r["protein_cluster"]
                     for r in csv.DictReader(fh, delimiter="\t")}
            missing = [c for c in cell_ids if c not in m]
            if missing:
                raise SystemExit(f"protein_clusters.tsv covers {len(m)} cells but "
                                 f"{len(missing)} obs cells are absent, e.g. {missing[:3]}")
            add_obs_column(f, "protein_cluster", [m[c] for c in cell_ids])
            print(f"merged protein_cluster ({len(set(m.values()))} clusters)")

    # 3. truth files
    with open(os.path.join(OUT, "pbmc.clusters_truth.tsv"), "w") as fh:
        fh.write("cell_id\ttruths\n")
        for c, l in zip(cell_ids, labels):
            fh.write(f"{c}\t{l}\n")
    n = len(set(labels))
    with open(os.path.join(OUT, "pbmc.clusters_truth_num.txt"), "w") as fh:
        fh.write(f"{n}\n")

    # 4. properties
    with open(os.path.join(OUT, "pbmc_properties.yaml"), "w") as fh:
        for k, v in PROPS.items():
            fh.write(f"{k}: {v}\n")

    # 5. checksums of the served files, in `sha256sum -c` format so it ports
    #    straight into the omni-data module (verify with `sha256sum -c`).
    served = ["pbmc.h5ad", "pbmc.clusters_truth.tsv",
              "pbmc.clusters_truth_num.txt", "pbmc_properties.yaml"]
    with open(os.path.join(OUT, "pbmc.sha256"), "w") as fh:
        for name in served:
            fh.write(f"{sha256(os.path.join(OUT, name))}  {name}\n")

    print(f"built {OUT}: {len(cell_ids)} cells, {n} {LABEL} labels")


if __name__ == "__main__":
    main()
