#/ Start module
module MDistributions

using Distributions
using Distributions: @check_args

using Optim
using Random
using RootSolvers
using StatsBase
using SpecialFunctions

#~ Pareto distributions
include("paretoI.jl")      #!note: special case, acts as wrapper that can be called as `Pareto`
include("paretoIV.jl")     #!note: also includes Pareto II and III as special cases
include("generalizedpareto.jl")
include("doublepareto.jl")
include("temperedpareto.jl")

#~ Tweedie distributions
include("tweedie.jl")

end # module Distributions
#/ End module




