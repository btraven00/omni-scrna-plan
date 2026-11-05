#!/usr/bin/env bash
# CI plan = the canonical plan with DATA swapped for the tiny committed fixture
# (omni-data reading fixtures/DATA via file://). Everything else is selected by
# a filter at run time (filters/ci-until-pca.yaml), not rewritten here.
#
# Usage: ci/fixture-plan.sh [SRC_YAML] > benchmark_ci.yaml
set -euo pipefail
src="${1:-benchmark_conda.yaml}"
root="$(cd "$(dirname "$0")/.." && pwd)"
export FIXTURE_URI="file://${root}/fixtures/DATA"
export OMNI_DATA_COMMIT="8df4cc8"   # omni-data main: --output id=PATH (api 0.8)

# Exactly one DATA stage, or the swap is a silent no-op and CI runs the real
# 817 MB dataset while calling itself a fixture run (the 2h45m incident).
if [[ "$(yq '[.stages[] | select(.id == "DATA")] | length' "$src")" != "1" ]]; then
  echo "fixture-plan: expected exactly one stage with id 'DATA' in $src" >&2
  exit 1
fi

# The fixture binds every DATA lineage label to "yes": unbound labels default
# to the module id, and FEAT/PCA/INTG8 require runs_*: "yes", so without this
# they are pruned and CI silently tests DATA..NORM only.
plan="$(yq '
  (.stages[] | select(.id == "DATA").modules) = [{
    "id": "be1",
    "name": "tiny be1 fixture (CI, file://)",
    "provides": ((.stages[] | select(.id == "DATA").provides) | map({"key": ., "value": "yes"}) | from_entries),
    "software_environment": "omni_data",
    "repository": {"url": "https://github.com/btraven00/omni-data", "commit": strenv(OMNI_DATA_COMMIT)},
    "parameters": [{"uri": strenv(FIXTURE_URI)}]
  }]
  | (.stages[] | select(.id == "DATA").modules[0].provides[]) style="double"  # "yes", not YAML 1.1 true
  | .software_environments.omni_data = {"name": "omni-data fixture loader", "conda": "envs/omni-data.yml"}
' "$src")"

if [[ "$(yq '.stages[] | select(.id == "DATA").modules[0].repository.url' <<<"$plan")" != *omni-data* ]]; then
  echo "fixture-plan: DATA swap did not take" >&2
  exit 1
fi
printf '%s\n' "$plan"
