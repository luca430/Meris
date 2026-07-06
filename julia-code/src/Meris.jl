module Meris

## GLOBAL CONSTANTS
const minsamples::Int = 30
const mincomponents::Int = 100

## DIRECTORIES
_with_trailing_separator(path) = joinpath(normpath(path), "")

const DATADIR = _with_trailing_separator(get(ENV, "MERIS_DATADIR", joinpath(@__DIR__, "..", "data")))
const FIGDIR = _with_trailing_separator(get(ENV, "MERIS_FIGDIR", joinpath(@__DIR__, "..", "figures")))

const ARXIVDIR = DATADIR * "datasets/arxiv/"
const TREEDIR = DATADIR * "datasets/bci.tree/"
const BIOTIMEDIR = DATADIR * "datasets/biotime/"
const BRIGHTKITEDIR = DATADIR * "datasets/brightkite/"
const EMAILDIR = DATADIR * "datasets/emails/"
const FINANCEDIR = DATADIR * "datasets/finance/"
const GAIADIR = DATADIR * "datasets/gaia/"
const GOWALLADIR = DATADIR * "datasets/gowalla/"
const GTEXDIR = DATADIR * "datasets/gtex/"
const GBIFDIR = DATADIR * "datasets/gbif/"
const GUTENBERGDIR = DATADIR * "datasets/gutenberg/"
const LEGODIR = DATADIR * "datasets/lego/"
const OTUDIR = DATADIR * "datasets/otu/"
const RFCDIR = DATADIR * "datasets/rfc/"
const TARADIR = DATADIR * "datasets/taraocean/"

## SUBMODULES
#~ CLI argument parser
include("args/argparse.jl")

#~ Data handlers, loaders, and samplers
include("dataframes/datatools.jl")

include("dataframes/arxiv/arxivloader.jl")
include("dataframes/bci.tree/bcitreeloader.jl")
include("dataframes/biotime/biotimeloader.jl")
include("dataframes/finance/financeloader.jl")
include("dataframes/gowalla/gowallaloader.jl")
include("dataframes/GTEx/gtexloader.jl")
include("dataframes/gutenberg/gutenbergloader.jl")
include("dataframes/lego/legoloader.jl")
include("dataframes/otu/otuloader.jl")
include("dataframes/rfc/rfcloader.jl")

#~ Processes
include("processes/dirichlet.jl")
include("processes/pitman-yor.jl")

#~ Utilities
include("distributions/distributions.jl")

include("dataframes/afd.jl")
include("dataframes/taylor.jl")

include("distributions/lr_distributions.jl")

#~ Fitting
include("fits/candidate-distributions.jl")
include("fits/goodness-of-fit.jl")
include("fits/straight-line.jl")
include("fits/heapsmodel.jl")

# include("fits/power-law.jl")
# include("fits/double-power-law.jl")
# include("fits/mle.jl")
# include("fits/heapsmodel.jl")

end # module Meris
