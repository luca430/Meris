#= Taylor's sampler =#
#/ Start module
module TaylorSampler

#/ Packages
using Distributions
using StatsBase
using Random
using SparseArrays

import Meris.MDistributions as MDist

#################
### FUNCTIONS ###

########################
### HELPER FUNCTIONS ###
function observe(
    N::Int, K::Int;
    γ::Float64=0.5,
    φ::Float64=1e5,
    ε::Float64=1e0,
    S::Int=2^14,
    rng=Random.Xoshiro(42*N)
)
    #~ Generate random propensities
    TPareto = MDist.TemperedPareto(γ, 1/φ, ε)
    w = MDist.rand(rng, TPareto, S)
    #~ Specify noise distribution
    Pξ = Distributions.Gamma()

    #~ Allocate
    x = K > 2N ? spzeros(Int, K, S) : zeros(Int, K, S)
    
    for k in 1:K
        ξ = rand(rng, Pξ, S)
        p = w .* ξ / sum(w .* ξ)
        n = StatsBase.sample(rng, 1:S, Weights(p), N; replace=true)
        cmap = StatsBase.countmap(n)
        for (label, counts) in cmap
            x[k, label] = counts
        end
    end
    μx, σx = taylor(x)
    return (; x=x, μx=μx, σx=σx)
end

########################
### HELPER FUNCTIONS ###
function taylor(x::Matrix{Int}; minoccupancy::Float64=0.1)
    K, S = size(x)
    #/ Compute mean and variance
    #~ Filter by occupancy
    occ = map(col -> count(!iszero, col), eachcol(x)) ./ K
    nzidx = findall(o -> o >= minoccupancy, occ)
    xnz = @view x[:,nzidx]
    μx = map(col -> mean(col), eachcol(xnz))
    σx = map(col -> var(col),  eachcol(xnz))
    return μx, σx
end

function taylor(x::SparseMatrixCSC; minoccupancy::Float64=0.1)
    K, S = size(x)
    #/ Compute mean and variance
    #~ Filter by occupancy
    occ = map(col -> count(!iszero, col), eachcol(x)) ./ K
    nzidx = findall(o -> o >= minoccupancy, occ)
    xnz = @view x[:,nzidx]
    μx = map(col -> mean(col), eachcol(xnz))
    σx = map(col -> var(col),  eachcol(xnz))
    return μx, σx
end

function ftaylor(x::Matrix{Int}, N::Int; samplingdepth::Float64 = 10.0)
    #/ Compute mean and variance
    μx = map(col -> mean(col), eachcol(x))
    σx = map(col -> var(col),  eachcol(x))
    fidxs = findall(μ -> μ > samplingdepth, μx)
    return μx[fidxs], σx[fidxs]
end

function ftaylor(x::SparseMatrixCSC, N::Int; samplingdepth::Float64 = 10.0)
    #/ Compute mean and variance
    μx = map(col -> mean(col), eachcol(x))
    σx = map(col -> var(col),  eachcol(x))
    fidxs = findall(μ -> μ > samplingdepth, μx)
    return μx[fidxs], σx[fidxs]
end

end # module TaylorSampler
#/ End module
