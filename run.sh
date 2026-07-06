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
#   RUN_INSTANTIATE=1
#   DATA_DIR=/data
#   RESULTS_DIR=/results
#   JULIA_FLAGS="--compiled-modules=no"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JULIA_PROJECT="$ROOT_DIR/julia-code"
if [[ -z "${RESULTS_DIR:-}" ]]; then
  if [[ -d "/results" || "$ROOT_DIR" == "/code" ]]; then
    RESULTS_DIR="/results"
  else
    RESULTS_DIR="$ROOT_DIR/results"
  fi
fi
REPRO_MODE="${REPRO_MODE:-bundled}"
FIGURE_SET="${FIGURE_SET:-main}"
CATEGORIES="${CATEGORIES:-linguistic,microbial,social,biology}"
JULIA_FLAGS="${JULIA_FLAGS:---compiled-modules=no}"
if [[ -n "${DATA_DIR:-}" ]]; then
  SOURCE_DATA_DIR="$DATA_DIR"
elif [[ -d "/data/macro" || -d "/data/fig3" ]]; then
  SOURCE_DATA_DIR="/data"
else
  SOURCE_DATA_DIR="$JULIA_PROJECT/data"
fi

WORK_DIR="$RESULTS_DIR/work/$REPRO_MODE"
WORK_DATA_DIR="$WORK_DIR/data"
PREDICTION_WORK_DIR="$WORK_DATA_DIR/macro/tl-prediction"
JULIA_RUNTIME_DEPOT="$RESULTS_DIR/julia_depot"
RUNTIME_PROJECT="$WORK_DIR/julia-project"

mkdir -p "$RESULTS_DIR/figures" "$RESULTS_DIR/tables" "$RESULTS_DIR/logs" "$PREDICTION_WORK_DIR" "$JULIA_RUNTIME_DEPOT" "$RUNTIME_PROJECT"

cp "$JULIA_PROJECT/Project.toml" "$RUNTIME_PROJECT/Project.toml"
ln -sfn "$JULIA_PROJECT/src" "$RUNTIME_PROJECT/src"

link_data_dir() {
  local relative_path="$1"
  local source_path="$SOURCE_DATA_DIR/$relative_path"
  local target_path="$WORK_DATA_DIR/$relative_path"

  if [[ -e "$source_path" && ! -e "$target_path" ]]; then
    mkdir -p "$(dirname "$target_path")"
    ln -s "$source_path" "$target_path"
  fi
}

link_prediction_file() {
  local category="$1"
  local source_path="$SOURCE_DATA_DIR/macro/tl-prediction/$category.jld2"
  local target_path="$PREDICTION_WORK_DIR/$category.jld2"

  if [[ ! -e "$source_path" ]]; then
    echo "Missing TL-prediction input: $source_path" >&2
    exit 4
  fi

  ln -sf "$source_path" "$target_path"
}

link_data_dir "macro/taylor"
link_data_dir "macro/sad"
link_data_dir "fig3"
link_data_dir "goodness-of-fit"
if [[ -d "$SOURCE_DATA_DIR/downsampled" ]]; then
  link_data_dir "downsampled"
fi

IFS=',' read -r -a CATEGORY_LIST <<< "$CATEGORIES"
if [[ "$REPRO_MODE" == "bundled" ]]; then
  for category in "${CATEGORY_LIST[@]}"; do
    link_prediction_file "$category"
  done
fi

export MERIS_DATADIR="$WORK_DATA_DIR"
export MERIS_FIGDIR="$JULIA_PROJECT/figures"
EXISTING_JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$HOME/.julia}"
export JULIA_DEPOT_PATH="$JULIA_RUNTIME_DEPOT:$EXISTING_JULIA_DEPOT_PATH:"

echo "Repository: $ROOT_DIR"
echo "Data:       $SOURCE_DATA_DIR"
echo "Results:    $RESULTS_DIR"
echo "Mode:       $REPRO_MODE"
echo "Figures:    $FIGURE_SET"
echo "Categories: $CATEGORIES"

if [[ "${RUN_INSTANTIATE:-0}" == "1" ]]; then
  echo "Instantiating Julia environment..."
  julia $JULIA_FLAGS --project="$RUNTIME_PROJECT" -e 'using Pkg; Pkg.instantiate()' \
    2>&1 | tee "$RESULTS_DIR/logs/instantiate.log"
fi

if ! julia $JULIA_FLAGS --project="$RUNTIME_PROJECT" -e 'using CSV, DataFrames, JLD2, LsqFit, CairoMakie' >/dev/null 2>&1; then
  echo "Julia dependencies are missing; instantiating project into the writable results depot..."
  julia --project="$RUNTIME_PROJECT" -e 'using Pkg; Pkg.instantiate()' \
    2>&1 | tee "$RESULTS_DIR/logs/instantiate-missing-dependencies.log"

  if ! julia $JULIA_FLAGS --project="$RUNTIME_PROJECT" -e 'using CSV, DataFrames, JLD2, LsqFit, CairoMakie' >/dev/null 2>&1; then
    echo "Julia dependencies are still unavailable after Pkg.instantiate()." >&2
    echo "See $RESULTS_DIR/logs/instantiate-missing-dependencies.log for details." >&2
    exit 10
  fi
fi

if [[ "$REPRO_MODE" == "regenerate" ]]; then
  if [[ ! -d "$WORK_DATA_DIR/downsampled" ]]; then
    echo "Missing optional grouped downsampled data: $SOURCE_DATA_DIR/downsampled" >&2
    echo "Use REPRO_MODE=bundled, or include data/downsampled to regenerate TL-prediction intermediates." >&2
    exit 3
  fi

  echo "Regenerating TL-prediction bins..."
  julia $JULIA_FLAGS --project="$RUNTIME_PROJECT" \
    "$JULIA_PROJECT/scripts/cli-scripts/macropatterns/tl-prediction.jl" \
    --categories="$CATEGORIES" \
    --downsampled-dir="$WORK_DATA_DIR/downsampled" \
    --result-dir="$PREDICTION_WORK_DIR" \
    2>&1 | tee "$RESULTS_DIR/logs/tl-prediction.log"
elif [[ "$REPRO_MODE" != "bundled" ]]; then
  echo "Unknown REPRO_MODE: $REPRO_MODE" >&2
  exit 2
fi

echo "Fitting power and quadratic Taylor-law candidates and computing AIC..."
julia $JULIA_FLAGS --project="$RUNTIME_PROJECT" \
  "$JULIA_PROJECT/scripts/cli-scripts/macropatterns/tl-prediction-aic.jl" \
  --categories="$CATEGORIES" \
  --prediction-dir="$PREDICTION_WORK_DIR" \
  2>&1 | tee "$RESULTS_DIR/logs/tl-prediction-aic.log"

cp "$PREDICTION_WORK_DIR/tl-prediction-aic.csv" \
  "$RESULTS_DIR/tables/"
cp "$PREDICTION_WORK_DIR/tl-prediction-aic-class-summary.csv" \
  "$RESULTS_DIR/tables/"

echo "Rendering figures..."
julia $JULIA_FLAGS --project="$RUNTIME_PROJECT" \
  "$ROOT_DIR/codeocean/reproduce_figures.jl" \
  --figure-set="$FIGURE_SET" \
  --results-dir="$RESULTS_DIR" \
  2>&1 | tee "$RESULTS_DIR/logs/reproduce-figures.log"

echo "Done. Reviewer outputs are in $RESULTS_DIR."
