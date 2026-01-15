#= Taylor's sampler =#
#/ Start module
module TaylorSampler

#/ Packages
using Distributions
using StatsBase
using Random
using SparseArrays
using JLD2

import Meris.MDistributions as MDist
import Meris.DATADIR as DATADIR
TLDIR = DATADIR * "taylor/synthetic/"

#################
### FUNCTIONS ###
"""
    observe

Observe system by multinomial sampling out of S components, each with prob. `p[i] = θ[i]/Θ`, where
`Θ = sum(θ)`, and `θ ~ TemperedPareto(γ,1/φ)` with `φ` the characteristic scale of the system.
, when `N≫φ` 
"""
function observe(
    N::Int, K::Int;
    γ::Float64=0.5,
    φ::Float64=1e5,
    ε::Float64=1e0,
    S::Int=10^5,
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
        # n = StatsBase.sample(rng, 1:S, Weights(p), N; replace=true)
        x[k,:] .+= rand.(rng, Poisson.(N.*p))
        
        # cmap = StatsBase.countmap(n)
        # for (label, counts) in cmap
        #     x[k, label] = counts > 1 ? counts : 0    #~ Measurement error, exclude hapax legomena
        # end
    end
    μx, σx = taylor(x)
    return (; x=x, μx=μx, σx=σx, params=(; N=N, K=K, S=S, φ=φ, γ=γ, ε=ε))
end

########################
### HELPER FUNCTIONS ###
function taylor(x::Matrix{Int}; minoccupancy::Float64=0., scale=false)
    K, S = size(x)
    #/ Compute mean and variance
    #~ Filter by occupancy
    occ = map(col -> count(!iszero, col), eachcol(x)) ./ K
    nzidx = findall(o -> o >= minoccupancy, occ)
    xnz = scale ? (@view x[:,nzidx]) : (@view x[:,:])
    μx = map(col -> mean(col), eachcol(xnz))
    σx = map(col -> var(col),  eachcol(xnz))
    return μx, σx
end

function taylor(x::SparseMatrixCSC; minoccupancy::Float64=0.01)
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

#############################
### DATA HELPER FUNCTIONS ###
function save_taylor(result, filename; DIR=TLDIR)
    mkpath(DIR)
    jldsave(DIR*filename, n=result.x, mean=result.μx, var=result.σx, params=result.params)
end

end # module TaylorSampler
#/ End module
