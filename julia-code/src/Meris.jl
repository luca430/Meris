module Meris

## DIRECTORIES
const DATADIR = normpath(joinpath(@__DIR__, "..", "data/"))
const BOOKDIR = DATADIR * "datasets/books/"
const LEGODIR = DATADIR * "datasets/lego/"
const OTUDIR = DATADIR * "datasets/otu/"
const ARXIVDIR = DATADIR * "datasets/arxiv/"
const CORPUSDIR = DATADIR * "datasets/corpus/"
const RFCDIR = DATADIR * "datasets/rfc/"
const TREEDIR = DATADIR * "datasets/bci.tree/"
const PATENTDIR = DATADIR * "datasets/patents/"
const GAIADIR = DATADIR * "datasets/gaia/"

## SUBMODULES
#~ Data handlers, loaders, and samplers 
include("dataframes/books/booksampler.jl")
include("dataframes/books/wordsampler.jl")
include("dataframes/lego/legosampler.jl")
include("dataframes/otu/otusampler.jl")
include("dataframes/arxiv/arxivsampler.jl")
include("dataframes/rfc/rfcsampler.jl")
include("dataframes/bci.tree/bctreesampler.jl")
include("dataframes/gaia/gaiasampler.jl")
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
# include("dataframes/taylor.jl")
include("dataframes/datatools.jl")
include("distributions/pareto.jl")
include("distributions/mle.jl")
include("distributions/lr_distributions.jl")

end # module Meris
