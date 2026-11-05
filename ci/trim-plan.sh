#!/usr/bin/env bash
# Derive a trimmed CI plan from the canonical benchmark_conda.yaml.
#
# The canonical plan stays the single source of truth; this computes a smaller
# view for CI smoke tests by:
#   1. swapping the DATA stage's download module for the tiny in-repo fixture
#      (omni-data reading fixtures/DATA via the file:// handler),
#   2. registering the omni_data software environment,
#   3. keeping only the first parameter combo-group per module, and
#   4. collapsing every remaining list-valued parameter to its first element
#      (so each module runs exactly once instead of a full sweep).
#
# Usage: ci/trim-plan.sh [SRC_YAML] > benchmark_ci.yaml
set -euo pipefail

src="${1:-benchmark_conda.yaml}"
root="$(cd "$(dirname "$0")/.." && pwd)"
fixture_uri="file://${root}/fixtures/DATA"

# Pinned omni-data commit (file:// fixture loader).
omni_data_commit="eb4d368b5e83b26389e6047dd8bad6badede78f8"

fixture_module="$(
  fixture_uri="$fixture_uri" omni_data_commit="$omni_data_commit" yq -n '
    [{
      "id": "be1",
      "name": "tiny be1 fixture (CI, non-reproducible file://)",
      "software_environment": "omni_data",
      "repository": {
        "url": "https://github.com/btraven00/omni-data",
        "commit": strenv(omni_data_commit)
      },
      "parameters": [{"uri": strenv(fixture_uri)}]
    }]'
)"

# The data download lives in the stage with id "DATA" (module id "data").
data_stage="DATA"

# Fail loudly if the stage we swap into doesn't exist — a missing match makes
# the assignment a silent no-op, so CI would run the FULL dataset while calling
# itself a "fixture run" (see the 2h45m / 1.5 GB-artifact incident).
if [[ "$(DS="$data_stage" yq '[.stages[] | select(.id == strenv(DS))] | length' "$src")" != "1" ]]; then
  echo "trim-plan: expected exactly one stage with id '$data_stage' in $src" >&2
  exit 1
fi

# The in_seurat RPCA module needs k_anchor=20 (not the collapsed-to-first 5) on
# the fixture: cross-type sample pairs (PBMC vs cell lines) yield few anchors, and
# at k_anchor=5 a pair drops below Seurat's k.weight and integration errors. 20 is
# an author-declared option ([5,20]) and clears it. Verified on real fixture inputs.
trimmed="$(
  FIXTURE_MODULE="$fixture_module" DS="$data_stage" yq '
    (.stages[] | select(.id == strenv(DS)).modules) = (strenv(FIXTURE_MODULE) | from_yaml)
    | .software_environments.omni_data = {"name": "omni-data fixture loader", "conda": "envs/omni-data.yml"}
    | (.stages[].modules[] | select(has("parameters")) | .parameters) |= [.[0]]
    | (.stages[].modules[] | select(has("parameters")) | .parameters[][] | select(tag == "!!seq")) |= .[0]
    | (.stages[] | select(.id == "INTG8").modules[] | select(.id == "in_seurat").parameters[0].k_anchor) = 20
  ' "$src"
)"

# Verify the swap actually landed: the DATA stage must now point at the fixture
# loader, not the real downloader. Guards against a future stage-id rename
# silently reintroducing the full-dataset run.
got_url="$(DS="$data_stage" yq '.stages[] | select(.id == strenv(DS)).modules[0].repository.url' <<<"$trimmed")"
if [[ "$got_url" != *omni-data* ]]; then
  echo "trim-plan: fixture swap did not take — DATA module repo is '$got_url', expected omni-data" >&2
  exit 1
fi

printf '%s\n' "$trimmed"

# NOTE: the canonical plan is api_version 0.4.0, which matches the pinned module
# commits (omnibenchmark passes --name=<dataset>, and modules build output
# filenames as "{name}_<suffix>"). Under api 0.5 --name becomes the module id and
# these modules would misname their outputs. Migrate the modules before bumping.
