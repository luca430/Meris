module Meris

## DIRECTORIES
const DATADIR   = normpath(joinpath(@__DIR__, "..", "data/"))
const BOOKDIR   = DATADIR * "datasets/books/"
const LEGODIR   = DATADIR * "datasets/lego/"
const OTUDIR    = DATADIR * "datasets/otu/"
const ARXIVDIR  = DATADIR * "datasets/arxiv/"
const CORPUSDIR = DATADIR * "datasets/corpus/"
const RFCDIR    = DATADIR * "datasets/rfc/"
const TREEDIR   = DATADIR * "datasets/bci.tree/"
const PATENTDIR = DATADIR * "datasets/patents/"
const WIKIDIR   = DATADIR * "datasets/wikitext-2/"

## SUBMODULES
#~ Data handlers, loaders, and samplers 
include("dataframes/books/booksampler.jl")
include("dataframes/books/wordsampler.jl")
include("dataframes/wikitext/wikitextsampler.jl")
include("dataframes/lego/legosampler.jl")
include("dataframes/otu/otusampler.jl")
include("dataframes/arxiv/arxivsampler.jl")
include("dataframes/rfc/rfcsampler.jl")
include("dataframes/bci.tree/bctreesampler.jl")
#~ Processes
include("processes/dirichlet.jl")
include("processes/pitman-yor.jl")
#~ Fitting
include("fits/goodness-of-fit.jl")
include("fits/straight-line.jl")
include("fits/power-law.jl")
include("fits/mle.jl")
#~ Utilities
include("dataframes/afd.jl")
include("dataframes/taylor.jl")
include("distributions/pareto.jl")
include("distributions/tweedie.jl")
end # module Meris
