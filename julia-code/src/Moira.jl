module Moira

const DATADIR = normpath(joinpath(@__DIR__, "..", "data/"))
const BOOKDIR = DATADIR * "datasets/books/"
const LEGODIR = DATADIR * "datasets/lego/"

## SUBMODULES
#~ Data samplers 
include("dataframes/books/booksampler.jl")
#~ Processes
include("processes/dirichlet.jl")
include("processes/pitman-yor.jl")


end # module Moira
