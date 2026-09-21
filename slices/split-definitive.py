#!/usr/bin/env python3
"""Generate the two-machine halves of definitive.yaml.

    python3 slices/split-definitive.py     # run from the repo root

The CPU half drops the GPU arms and gains pc-scrapper; the GPU half keeps
only the rapids arms plus the upstream needed to feed them. Everything the
two share is emitted byte-identical, so the param hashes agree and the two
out/ trees merge by path. Local-path module pins are rewritten to the
btraven00 forks (a remote host cannot resolve /home/b/...), and the -prof
entrypoints to plain ones (they invoke denet, a cargo binary on neither host).
"""
import re, pathlib

src = pathlib.Path("slices/definitive.yaml").read_text().split("\n")

MOD = re.compile(r"^      - id: (\S+)\s*$")
STAGE = re.compile(r"^  - id: (\S+)\s*$")

def blocks(lines):
    """(name, start, end) per module block; start absorbs the comment lines above it."""
    starts = [(i, m.group(1)) for i, l in enumerate(lines) if (m := MOD.match(l))]
    out = []
    for n, (i, name) in enumerate(starts):
        s = i
        while s > 0 and (lines[s - 1].startswith("      #")):
            s -= 1
        # end: next module start (with its comments) or the next stage/separator
        e = len(lines)
        for j in range(i + 1, len(lines)):
            if MOD.match(lines[j]) or STAGE.match(lines[j]) or lines[j].startswith("  # ---"):
                e = j
                break
        while e > i and lines[e - 1].strip() == "":
            e -= 1
        out.append((name, s, e))
    # re-resolve starts so a block's absorbed comments do not overlap the one above
    fixed = []
    for k, (name, s, e) in enumerate(out):
        if k and s < out[k - 1][2]:
            s = out[k - 1][2]
        fixed.append((name, s, e))
    return fixed

def drop(lines, names):
    keep, bl = [], blocks(lines)
    doomed = set()
    for name, s, e in bl:
        if name in names:
            doomed.update(range(s, e))
            # swallow the blank line that followed the block
            if e < len(lines) and lines[e].strip() == "":
                doomed.add(e)
    return [l for i, l in enumerate(lines) if i not in doomed]

def insert_after(lines, name, text):
    for n, s, e in blocks(lines):
        if n == name:
            return lines[:e] + [""] + text.split("\n") + lines[e:]
    raise SystemExit(f"anchor {name} not found")

PC_SCRAPPER = '''      # Third PCA arm, and the only one that threads: libscran's irlba behind
      # scrapper::runPca, with --num_threads wired through. scanpy's randomized
      # arm needs dense input to have a seed axis at all (it silently coerces to
      # arpack on sparse), so this is the sparse-native seed axis, and the one
      # that can actually use the 64 cores on the CPU host.
      - id: pc-scrapper
        name: "scrapper PCA (libscran irlba, threaded)"
        software_environment: "scrapper"
        repository:
          url: https://github.com/btraven00/scrapper
          commit: 4f14524
          entrypoint: pca
        parameters:
          - solver: irlba
            backend: memory        # matched to scanpy/rapids: no DelayedArray
            dense: "false"
            num_threads: 8
            n_components: 50
            random_seed: *pca_seeds'''

def retitle(lines, new_id, blurb):
    out = list(lines)
    for i, l in enumerate(out[:12]):
        if l.startswith("name: "):
            out[i] = f"name: {new_id}"
        elif l.startswith("id: "):
            out[i] = f"id: {new_id.replace('-', '_')}"
    return [blurb] + out

def unprof(lines):
    return [l.replace("entrypoint: pca-prof", "entrypoint: pca")
             .replace("entrypoint: cluster-prof", "entrypoint: cluster") for l in lines]


STAGE_RE = __import__("re").compile(r"^  - id: (\S+)\s*$")

def drop_stage(lines, name):
    """Remove a whole stage, plus the separator comment block above it."""
    at = next(i for i, l in enumerate(lines) if STAGE_RE.match(l) and STAGE_RE.match(l).group(1) == name)
    start = at
    while start > 0 and (lines[start - 1].startswith("  #") or lines[start - 1].strip() == ""):
        start -= 1
    end = len(lines)
    # from `at`, not `start`: start has been walked back over the comment block,
    # so searching from start+1 would find this very stage and delete nothing.
    for j in range(at + 1, len(lines)):
        if STAGE_RE.match(lines[j]):
            end = j
            break
    while end > start and lines[end - 1].strip() == "":
        end -= 1
    keep = lines[:start] + lines[end:]
    return keep

REPIN = [
    ("          url: /home/b/lab/omni-scrna/split-stages-plan-modules/sc3s\n          commit: 502f6da",
     "          url: https://github.com/btraven00/sc3s\n          commit: e6e6a85"),
    ("          url: /home/b/lab/omni-scrna/split-stages-plan-modules/metrics\n          commit: bc998fb",
     "          url: https://github.com/btraven00/omni-scrna-metrics\n          commit: 9aedfb5"),
]

def repin(lines):
    t = "\n".join(lines)
    for a, b in REPIN:
        t = t.replace(a, b)
    assert "url: /home/b" not in t, "a local-path pin survived"
    return t.split("\n")

HDR = ("# GENERATED from definitive.yaml -- edit that, then regenerate.\n"
       "# Split for a two-machine run; the shared DATA/FILT/NORM/FEAT/NNG stanzas\n"
       "# are byte-identical to the GPU half so the param hashes match and the two\n"
       "# out/ trees merge by path.\n"
       "# The -prof entrypoints are swapped for plain ones: they invoke denet, which\n"
       "# is a cargo binary present on neither host.\n")

HDR_GPU = ("# GENERATED from definitive.yaml -- edit that, then regenerate.\n"
           "# The rapids half of the two-machine split: only the GPU arms, plus the\n"
           "# shared upstream needed to feed them (byte-identical to the CPU half).\n"
           "# CLUST-E/sc3s is CPU work and lives in the other slice.\n"
           "# The -prof entrypoints are swapped for plain ones: they invoke denet,\n"
           "# which is a cargo binary present on neither host.\n")

cpu = repin(unprof(src))
cpu = drop(cpu, {"pc-rapids", "cl-rapids"})
cpu = insert_after(cpu, "pc-scanpy", PC_SCRAPPER)
cpu = retitle(cpu, "definitive-cpu", HDR)
pathlib.Path("slices/definitive-cpu.yaml").write_text("\n".join(cpu))

gpu = drop_stage(repin(unprof(src)), 'CLUST-E')
gpu = drop(gpu, {"pc-scanpy", "cl-scanpy", "cl-scrapper", "cl-seurat", "cl-sc3s"})
gpu = retitle(gpu, "definitive-gpu", HDR_GPU)
pathlib.Path("slices/definitive-gpu.yaml").write_text("\n".join(gpu))
print("wrote both")

# --- per-dataset CPU halves -------------------------------------------------
# One dataset at a time, so a run can be placed on a machine that suits its
# size. The DATA stanzas are untouched, so every node hashes the same as it
# would in the combined slice and the out/ trees still merge by path.
DATASETS = ["d-integration", "tm-facs", "tm-droplet", "aifi-d1"]

for ds in DATASETS:
    one = drop(list(cpu), {d for d in DATASETS if d != ds})
    one = [l.replace("name: definitive-cpu", f"name: definitive-cpu-{ds}")
            .replace("id: definitive_cpu", f"id: definitive_cpu_{ds.replace('-', '_')}")
           for l in one]
    pathlib.Path(f"slices/definitive-cpu-{ds}.yaml").write_text("\n".join(one))
print("wrote per-dataset slices:", ", ".join(DATASETS))

# --- rapids top-up ----------------------------------------------------------
# cl-rapids was severed from the seed sweep by the CPU/GPU split: it lived in
# the GPU half, which had no pc-scanpy or pc-scrapper, so it only ever attached
# to the single deterministic pc-rapids embedding -- 1 PCA seed instead of 20,
# no crossed design, and below pairwise.py's --min-runs.
#
# This slice re-attaches it: the CPU half's PCA stanzas verbatim (so the 40
# existing embeddings are reused, not recomputed) with cl-rapids as the only
# clusterer. CLUST-E and CLUSTBOUND are dropped -- already computed, and
# leaving them out keeps the DAG small.
CL_RAPIDS = [l for l in blocks(gpu) if l[0] == "cl-rapids"]
_n, _s, _e = CL_RAPIDS[0]
rapids_block = gpu[_s:_e]

top = drop(list(cpu), {"cl-scanpy", "cl-scrapper", "cl-seurat"})
for st in ("CLUST-E", "CLUSTBOUND"):
    top = drop_stage(top, st)
# put cl-rapids into the now-empty CLUST stage
out, i = [], 0
while i < len(top):
    out.append(top[i])
    if top[i].strip() == "modules:" and any(
            l.strip() == "- id: CLUST" for l in top[max(0, i - 12):i]):
        out.extend(rapids_block)
    i += 1
out = [l.replace("name: definitive-cpu", "name: definitive-rapids-topup")
        .replace("id: definitive_cpu", "id: definitive_rapids_topup") for l in out]
pathlib.Path("slices/definitive-rapids-topup.yaml").write_text("\n".join(out))
print("wrote slices/definitive-rapids-topup.yaml")
