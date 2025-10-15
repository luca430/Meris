module Meris

## DIRECTORIES
const DATADIR   = normpath(joinpath(@__DIR__, "..", "data/"))
const BOOKDIR   = DATADIR * "datasets/books/"
const LEGODIR   = DATADIR * "datasets/lego/"
const OTUDIR    = DATADIR * "datasets/otu/"
const ARXIVDIR  = DATADIR * "datasets/arxiv/"
const CORPUSDIR = DATADIR * "datasets/corpus/"
const RFCDIR    = DATADIR * "datasets/rfc/"

## SUBMODULES
#~ Data samplers 
include("dataframes/books/booksampler.jl")
include("dataframes/books/wordsampler.jl")
include("dataframes/lego/legosampler.jl")
include("dataframes/otu/otusampler.jl")
include("dataframes/arxiv/arxivsampler.jl")
#~ Processes
include("processes/dirichlet.jl")
include("processes/pitman-yor.jl")
#~ Utilities
include("dataframes/afd.jl")
include("dataframes/taylor.jl")
include("distributions/mle.jl")

end # module Meris
