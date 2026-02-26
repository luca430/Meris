# compositional

Repository for studying the emergence of statistical laws in complex component systems.

The main Julia package, `Meris`, provides tools for:

- dataset loading and preprocessing,
- analysis of Taylor's law,
- analysis and fitting of heavy-tailed component abundance distributions (CADs),
- analysis and fitting of Heaps' law,
- candidate-model comparison and goodness-of-fit workflows.

The repository also includes Python and Bash scripts to download the datasets used in the analyses.

## Project structure

- `julia-code/src/Meris.jl`: package entrypoint and module wiring.
- `julia-code/src/dataframes/`: dataset loaders and dataframe utilities.
- `julia-code/src/fits/`: fitting routines, candidate distributions, goodness-of-fit.
- `julia-code/src/distributions/`: custom distribution implementations.
- `julia-code/src/processes/`: stochastic process models.
- `julia-code/scripts/cli-scripts/`: CLI analysis pipelines.
- `julia-code/data/datasets/`: dataset download/parsing scripts.
- `python-notebooks/`: exploratory jupyter notebooks (using julia code).

## Setup (Julia)

```bash
cd julia-code
julia --project -e "using Pkg; Pkg.instantiate()"
```

Then in Julia:

```julia
using Meris
```

## Run analysis scripts

Most analysis entrypoints are in:

- `julia-code/scripts/cli-scripts/scaling-laws/`
- `julia-code/scripts/cli-scripts/goodness-of-fit/`
- `julia-code/scripts/cli-scripts/macro-laws/`

Example:

```bash
cd julia-code
julia --project scripts/cli-scripts/scaling-laws/heaps/rfc-heaps.jl
```

## Dataset download scripts

Dataset helpers are in `julia-code/data/datasets/`, with ready-to-run examples in many `command.txt` files.

Run each downloader from its own dataset folder.
