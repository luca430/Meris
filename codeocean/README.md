# CodeOcean Reproduction Workflow

Use the repository root `run.sh` as the CodeOcean run script.

Default command:

```bash
./run.sh
```

This default, `REPRO_MODE=bundled`, is the reviewer-friendly path. It uses the
checked-in intermediate analysis data under `julia-code/data/`, recomputes the
Taylor-law model comparison for the power and quadratic forms, writes AIC
tables, and renders the main figures.

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

The regenerate mode uses the bundled files in `julia-code/data/downsampled/`.
It does not download external raw datasets.

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

- Taylor-law summaries: `julia-code/scripts/module-scripts/macropatterns/taylor.jl`
- TL-prediction bins: `julia-code/scripts/cli-scripts/macropatterns/tl-prediction.jl`
- Power/quadratic AIC comparison: `julia-code/scripts/cli-scripts/macropatterns/tl-prediction-aic.jl`
- Goodness-of-fit candidate fits: `julia-code/scripts/cli-scripts/goodness-of-fit/`
- Figure scripts: `julia-code/plot/`
