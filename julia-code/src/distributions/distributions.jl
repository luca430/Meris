#/ Start module
module MDistributions

using Distributions
using Distributions: @check_args

using Optim
using Random
using StatsBase


#~ Include distributions
include("paretoI.jl")
include("doublepareto.jl")

end # module Distributions
#/ End module




