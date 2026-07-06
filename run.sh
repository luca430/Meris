#!/usr/bin/env bash
set -euo pipefail

# CodeOcean/reviewer entrypoint.
#
# Default mode uses the checked-in intermediate JLD2 files, recomputes the
# Taylor-law candidate AIC table, and redraws the main paper figures.
#
# Optional environment variables:
#   REPRO_MODE=bundled|regenerate
#   FIGURE_SET=main|taylor|all
#   CATEGORIES=linguistic,microbial,social,biology
#   SKIP_INSTANTIATE=1
#   RESULTS_DIR=/results

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JULIA_PROJECT="$ROOT_DIR/julia-code"
RESULTS_DIR="${RESULTS_DIR:-$ROOT_DIR/results}"
REPRO_MODE="${REPRO_MODE:-bundled}"
FIGURE_SET="${FIGURE_SET:-main}"
CATEGORIES="${CATEGORIES:-linguistic,microbial,social,biology}"

mkdir -p "$RESULTS_DIR/figures" "$RESULTS_DIR/tables" "$RESULTS_DIR/logs"

echo "Repository: $ROOT_DIR"
echo "Results:    $RESULTS_DIR"
echo "Mode:       $REPRO_MODE"
echo "Figures:    $FIGURE_SET"
echo "Categories: $CATEGORIES"

if [[ "${SKIP_INSTANTIATE:-0}" != "1" ]]; then
  echo "Instantiating Julia environment..."
  julia --project="$JULIA_PROJECT" -e 'using Pkg; Pkg.instantiate()' \
    2>&1 | tee "$RESULTS_DIR/logs/instantiate.log"
fi

if [[ "$REPRO_MODE" == "regenerate" ]]; then
  echo "Regenerating grouped downsampled datasets..."
  julia --project="$JULIA_PROJECT" \
    "$JULIA_PROJECT/scripts/cli-scripts/downsample-dataset-groups.jl" \
    --groups="$CATEGORIES" \
    2>&1 | tee "$RESULTS_DIR/logs/downsample-dataset-groups.log"

  echo "Regenerating TL-prediction bins..."
  julia --project="$JULIA_PROJECT" \
    "$JULIA_PROJECT/scripts/cli-scripts/macropatterns/tl-prediction.jl" \
    --categories="$CATEGORIES" \
    2>&1 | tee "$RESULTS_DIR/logs/tl-prediction.log"
elif [[ "$REPRO_MODE" != "bundled" ]]; then
  echo "Unknown REPRO_MODE: $REPRO_MODE" >&2
  exit 2
fi

echo "Fitting power and quadratic Taylor-law candidates and computing AIC..."
julia --project="$JULIA_PROJECT" \
  "$JULIA_PROJECT/scripts/cli-scripts/macropatterns/tl-prediction-aic.jl" \
  --categories="$CATEGORIES" \
  2>&1 | tee "$RESULTS_DIR/logs/tl-prediction-aic.log"

cp "$JULIA_PROJECT/data/macro/tl-prediction/tl-prediction-aic.csv" \
  "$RESULTS_DIR/tables/"
cp "$JULIA_PROJECT/data/macro/tl-prediction/tl-prediction-aic-class-summary.csv" \
  "$RESULTS_DIR/tables/"

echo "Rendering figures..."
julia --project="$JULIA_PROJECT" \
  "$ROOT_DIR/codeocean/reproduce_figures.jl" \
  --figure-set="$FIGURE_SET" \
  --results-dir="$RESULTS_DIR" \
  2>&1 | tee "$RESULTS_DIR/logs/reproduce-figures.log"

echo "Done. Reviewer outputs are in $RESULTS_DIR."
