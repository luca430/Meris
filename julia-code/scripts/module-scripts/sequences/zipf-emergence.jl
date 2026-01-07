#= Zipf's generator =#
#/ Start module
module ZipfGenerator

#/ Packages
using Distributions
using StatsBase
using Random

import Meris.MDistributions as MDist

#################
### FUNCTIONS ###
function observe(
    N::Int;
    γ::Float64=0.5,
    φ::Float64=1e5,
    ε::Float64=1e0,
    S::Int=2^12,
    rng=Random.Xoshiro(42*N)
)
    #~ Generate random propensities
    pdf = MDist.ParetoI(γ, ε)
    w = MDist.rand(rng, pdf, S)
    #~ Specify noise distribution
    Pξ = Distributions.Gamma()

    x = zeros(Int, S)

    #~ Allocate    
    # ξ = rand(rng, Pξ, S)
    # p = w .* ξ / sum(w .* ξ)
    p = w / sum(w)
    return w, p
    n = StatsBase.sample(rng, 1:S, Weights(p), N; replace=true)
    cmap = StatsBase.countmap(n)
    for (label, counts) in cmap
        x[label] = counts
    end

    Z = zipf(x)
    return (; p=p, n=Z.n)
end

########################
### HELPER FUNCTIONS ###
function zipf(x::Array{Int})
    #~ Compute sorted counts
    xnz = x[x .> 0]
    ns = xnz ./ sum(xnz)
    ranks = sortperm(ns, rev=true)
    return (; n=ns[ranks])
end

end # module ZipfGenerator
#/ End module
