#!/usr/bin/env bash
# Run the bundled fixture through scSHC. Override H5 / TSVs to point at real data.
set -euo pipefail
cd "$(dirname "$0")"

H5="${H5:-data/data.h5ad}"
OUT="${OUT:-out}"
CORES="${SCSHC_CORES:-8}"

# default: the bundled real cluster tsv; pass extra tsv paths as args to sweep many
TSVS=("$@")
[ ${#TSVS[@]} -eq 0 ] && TSVS=(data/*.tsv)

SCSHC_CORES="$CORES" Rscript scshc_eval.R "$H5" "$OUT" "${TSVS[@]}"
echo "--- summary ---"; cat "$OUT/summary.tsv"
