module Meris

## DIRECTORIES
const DATADIR = normpath(joinpath(@__DIR__, "..", "data/"))
const BOOKDIR = DATADIR * "datasets/books/"
const LEGODIR = DATADIR * "datasets/lego/"
const OTUDIR  = DATADIR * "datasets/otu/"

const CORPUSDIR = DATADIR * "datasets/corpus/"

## SUBMODULES
#~ Data samplers 
include("dataframes/books/booksampler.jl")
include("dataframes/lego/legosampler.jl")
include("dataframes/otu/otusampler.jl")
#~ Processes
include("processes/dirichlet.jl")
include("processes/pitman-yor.jl")

end # module Meris
