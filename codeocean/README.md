# CodeOcean Reproduction Workflow

Use the repository root `run.sh` as the CodeOcean run script.

Default command:

```bash
cd /code
RESULTS_DIR=/results ./run.sh
```

This default, `REPRO_MODE=bundled`, is the reviewer-friendly path. It uses the
prepared intermediate analysis data under `/data`, recomputes the
Taylor-law model comparison for the power and quadratic forms, writes AIC
tables, and renders the main figures.

The bundled mode does **not** require the heavy grouped downsampled data in
`/data/downsampled/`. It only needs the prepared intermediates in
`/data/macro/`, `/data/fig3/`, and `/data/goodness-of-fit/`.

Install Julia and the package dependencies with the CodeOcean environment
editor before runtime. The run script does not install packages unless
`RUN_INSTANTIATE=1` is set explicitly for a local/debug run.

If CodeOcean only offers Julia 1.12, upload `julia-code/Project.toml` but do
not upload the local `julia-code/Manifest.toml` generated with Julia 1.11.
The `environment/Project.toml` file mirrors the project dependencies so the
CodeOcean Docker build can instantiate them with Julia 1.12 before `/code` is
mounted. The run script then creates a writable runtime Julia project under
`/results/work/` and uses those installed dependencies.

In this CodeOcean interface, `environment/postInstall` runs before `/code` is
mounted, so it cannot access `julia-code/Project.toml`. The runtime script
therefore checks for missing Julia dependencies and instantiates the project
into `/results/julia_depot` if needed.

Outputs are written to:

- `results/tables/tl-prediction-aic.csv`
- `results/tables/tl-prediction-aic-class-summary.csv`
- `results/figures/fig2_A.pdf`
- `results/figures/fig2_B.pdf`
- `results/figures/fig3.pdf`
- `results/logs/`

## Reproduction Modes

Fast/default mode:

```bash
REPRO_MODE=bundled ./run.sh
```

Regenerate TL-prediction intermediates from the grouped downsampled datasets
before recomputing AIC and figures:

```bash
REPRO_MODE=regenerate ./run.sh
```

The regenerate mode requires the optional grouped downsampled files in
`/data/downsampled/`:

- `linguistic.jld2`
- `microbial.jld2`
- `social.jld2`
- `biology.jld2`

These files are large and may be omitted from a lightweight CodeOcean capsule.
If they are omitted, use `REPRO_MODE=bundled`. Regenerate mode does not
download external raw datasets; it rebuilds TL-prediction intermediates from
those grouped downsampled files.

## Figure Sets

Main manuscript figures:

```bash
FIGURE_SET=main ./run.sh
```

Taylor-law and TL-prediction supplementary figures:

```bash
FIGURE_SET=taylor ./run.sh
```

Both sets:

```bash
FIGURE_SET=all ./run.sh
```

## Dataset Downloaders

Raw dataset download helpers are under `julia-code/data/datasets/`. Several
folders contain a `command.txt` with the original command. Those scripts are
kept separate from the default CodeOcean run because they depend on external
services and can be slow or change outside the capsule.

## Main Analysis Entry Points

- TL-prediction bins: `julia-code/scripts/cli-scripts/macropatterns/tl-prediction.jl`
- Power/quadratic AIC comparison: `julia-code/scripts/cli-scripts/macropatterns/tl-prediction-aic.jl`
- Grouped downsampled data regeneration: `julia-code/scripts/cli-scripts/downsample-dataset-groups.jl`
- SAD helper used by Figure 3: `julia-code/scripts/module-scripts/macropatterns/SAD.jl`
- Figure scripts: `julia-code/plot/`
