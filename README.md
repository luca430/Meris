# Meris

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
- `julia-code/data/downsampled/`: grouped downsampled datasets used by downsampled plots.
- `julia-code/data/macro/`: checked-in intermediate results for Taylor's law, SAD, and TL-prediction plots.
- `julia-code/data/fig3/`: checked-in intermediate results for Figure 3.
- `julia-code/plot/`: figure-generation scripts.
- `julia-code/figures/`: generated figure PDFs and image assets.
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

Unless otherwise noted, the commands below are run from the repository root.

## CodeOcean / reviewer reproduction

The root `run.sh` script is the recommended CodeOcean entrypoint:

```bash
./run.sh
```

By default it uses the checked-in intermediate data under `julia-code/data/`,
recomputes the Taylor-law power-vs-quadratic AIC table, and redraws the main
figures into `results/`.

Useful variants:

```bash
REPRO_MODE=regenerate ./run.sh   # rebuild TL-prediction bins from grouped downsampled data
FIGURE_SET=taylor ./run.sh       # Taylor-law and TL-prediction supplementary figures
FIGURE_SET=all ./run.sh          # main + Taylor-law figure sets
```

See `codeocean/README.md` for the full capsule workflow and output list.

## Load or import data

Raw dataset download and parsing helpers live under `julia-code/data/datasets/`.
Several dataset folders include a `command.txt` with a ready-to-run example.
Run those commands from the dataset folder itself, for example:

```bash
cd julia-code/data/datasets/rfc
./download.sh
```

After raw files are present, use the Julia loaders exposed by `Meris`. They
return standardized `DataFrame`s with columns such as `class`, `sample_id`,
`component_id`, `counts`, and `nreads`.

```julia
using Meris

df_rfc = Meris.RFCLoader.load()
df_otu = Meris.OTULoader.load()
df_lego = Meris.LEGOLoader.load()
```

The main checked-in analysis data are JLD2 files under `julia-code/data/`.
For example, the Taylor-law and TL-prediction inputs used by Figure 2 are in:

- `julia-code/data/macro/taylor/`
- `julia-code/data/macro/tl-prediction/`

Figure 3 reads prepared SAD intermediates from:

- `julia-code/data/fig3/`

You can inspect those files directly:

```julia
using JLD2

taylor = load("julia-code/data/macro/taylor/linguistic.jld2")
prediction = load("julia-code/data/macro/tl-prediction/linguistic.jld2")
fig3_part = load("julia-code/data/fig3/linguistic/rfc.jld2")
```

## Downsampled data

Grouped downsampled datasets are stored in `julia-code/data/downsampled/`.
Each group file contains `ds_df` plus summary tables:

- `linguistic.jld2`: arXiv, Gutenberg, RFC
- `microbial.jld2`: OTU
- `social.jld2`: finance, Gowalla, LEGO
- `biology.jld2`: GTEx, BCI trees, BioTIME

Load a grouped downsampled dataset with:

```julia
using JLD2

data = load("julia-code/data/downsampled/linguistic.jld2")
df = data["ds_df"]
summary = data["summary"]
```

To regenerate the grouped downsampled files:

```bash
julia --project=julia-code julia-code/scripts/cli-scripts/downsample-dataset-groups.jl
```

Useful options:

```bash
julia --project=julia-code julia-code/scripts/cli-scripts/downsample-dataset-groups.jl --groups=linguistic,social --seed=123 --csv
```

To downsample each source dataset separately instead of by plotted group:

```bash
julia --project=julia-code julia-code/scripts/cli-scripts/downsample-datasets.jl
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

## Generate main figures

The repository includes generated PDFs in `julia-code/figures/`:

- `fig2_A.pdf`
- `fig2_B.pdf`
- `fig3.pdf`

To regenerate Figure 2A and Figure 2B from the checked-in Taylor-law and
TL-prediction intermediates:

```bash
julia --project=julia-code -e 'include("julia-code/plot/plot-fig2-A.jl"); Figure2A.plot()'
julia --project=julia-code -e 'include("julia-code/plot/plot-fig2-B.jl"); Figure2B.plot()'
```

To regenerate Figure 3 from the checked-in `julia-code/data/fig3/`
intermediates:

```bash
julia --project=julia-code -e 'include("julia-code/plot/plot-fig3.jl"); Figure3.plot()'
```

If the Figure 3 intermediates need to be rebuilt from raw loaders and
goodness-of-fit outputs, run:

```bash
julia --project=julia-code -e 'include("julia-code/plot/plot-fig3.jl"); Figure3.prepare(); Figure3.plot()'
```

`Figure3.prepare()` expects the standardized raw datasets to be available
through the `Meris.*Loader.load()` functions and candidate-fit files in
`julia-code/data/goodness-of-fit/`.
