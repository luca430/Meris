module Meris

## DIRECTORIES
const DATADIR = normpath(joinpath(@__DIR__, "..", "data/"))
const FIGDIR = normpath(joinpath(@__DIR__, "..", "figures/"))

const ARXIVDIR = DATADIR * "datasets/arxiv/"
const TREEDIR = DATADIR * "datasets/bci.tree/"
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
#~ Data handlers, loaders, and samplers
include("dataframes/arxiv/arxivloader.jl")
include("dataframes/bci.tree/bcitreesampler.jl")
include("dataframes/gaia/gaialoader.jl")
include("dataframes/gowalla/gowallaloader.jl")
include("dataframes/GTEx/gtexloader.jl")
# include("dataframes/gutenberg/gutenbergloader.jl")
include("dataframes/lego/legosampler.jl")
include("dataframes/otu/otuloader.jl")
include("dataframes/rfc/rfcsampler.jl")
# include("dataframes/tara/taraloader.jl")

#~ Processes
include("processes/dirichlet.jl")
include("processes/pitman-yor.jl")

#~ Fitting
include("fits/goodness-of-fit.jl")
include("fits/straight-line.jl")
include("fits/power-law.jl")
include("fits/double-power-law.jl")
include("fits/mle.jl")

#~ Utilities
include("dataframes/afd.jl")
include("dataframes/taylor.jl")
include("dataframes/datatools.jl")
include("distributions/distributions.jl")
include("distributions/lr_distributions.jl")
end # module Meris
