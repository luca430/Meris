#= Relatively simple module with function to illustrate that Taylor's law with
   exponent 2 emerges naturally from a Dirichlet prior with heavy-tailed pseudo-counts
=#
#/ Start module
module SaturatedTaylor

#/ Packages
using Distributions
using StatsBase
using Random

import Meris.MDistributions as MDist
import Meris.PitmanYor as PY

#################
### FUNCTIONS ###
"""Obtain random [ordered] pseudo-counts

!note: Note that ordering preserves the expected value of α[i], as it mimicks an atom distribution
@TODO: Show this
"""
function generate_pseudocounts(γ::Float64, φ::Float64; S::Int=2^14, ε=1.0, rng=Random.Xoshiro(42))
	  P = MDist.TemperedPareto(γ, 1/φ, ε)
    α = MDist.rand(rng, P, S)
    p = α ./ sum(α)
    return (α, p)
end

function generate_fixedpseudocounts(φ::Float64, γ::Float64; S::Int=2^14)
    αv = map(k -> k^(-γ), 1:S)
    Z = sum(αv)
    α = (φ / Z) .* αv
    return α
end

function generate_pseudocounts(S::Int; θ::Float64=1e0, α::Float64=0.25, rng=Random.Xoshiro(42))
    _, α = PY.generate_pseudocounts(S; θ=θ, α=α, rng=rng)
    p = α ./ sum(α)    
    return (α, p)
end

"""Obtain random observations of length N"""
function generate_observation(N::Int, p::Vector{Float64}; rng=Random.Xoshiro(42*N))
    S = length(p)
    n = StatsBase.sample(rng, 1:S, Weights(p), N; replace=true)
    return n
end

"""Observe"""
function observe(K::Int, N::Int; γ::Float64=0.8, φ::Float64=1e5, S::Int=2^10)
    x = zeros(Int, K, S)
    φv = Array{Float64}(undef, K)
    rng = Random.Xoshiro(31*K*N)
    
    Threads.@threads for k in 1:K
        #~ Generate system specific properties
        # α, p = generate_pseudocounts(γ, φ; S=S, rng=rng)
        α, p = generate_pseudocounts(S; α=0.2)
        φv[k] = sum(α)
        #~ Observe system by taking N samples
        n = generate_observation(N, p; rng=rng)
        #~ Count occurences
        cmap = StatsBase.countmap(n)
        for (label, counts) in cmap
            x[k, label] = counts
        end
    end
    return (; x=x, φ=φv)
end

function observe_Dirichlet(
    K::Int, N::Int;
    γ::Float64=1.5,
    μφ::Float64=1e1,
    σφ::Float64=1e0,
    S::Int=2^12
)
    x = zeros(Int, K, S)
    rng = Random.Xoshiro(42*K*N)
    #~ Generate pseudocounts
    # α = generate_fixedpseudocounts(φ, γ; S=S)
    # Dir = Distributions.Dirichlet(α)

    for k in 1:K
        φk = rand(rng, LogNormal(μφ, σφ))
        α = generate_fixedpseudocounts(φk, γ; S=S)
        p = rand(rng, Dirichlet(α))

        #~ Generate probabilities
        # p = rand(rng, Dir)
        #~ Observe system by taking N samples
        n = generate_observation(N, p; rng=rng)
        #~ Count occurences
        cmap = StatsBase.countmap(n)
        for (label, counts) in cmap
            x[k, label] = counts
        end
    end
    μx, σx = ftaylor(x, N)
    return (; x=x, μx=μx, σx=σx)
end

function taylor(x::Matrix{Int}; minoccupancy::Float64=1.0)
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
    K, S = size(x)
    #~ Normalize
    # xv = reduce(hcat, map(row -> row / sum(row), eachrow(x)))
    #/ Compute mean and variance
    μx = map(col -> mean(col), eachcol(x))
    σx = map(col -> var(col),  eachcol(x))
    fidxs = findall(μ -> N*μ > samplingdepth, μx)
    return μx[fidxs], σx[fidxs]
end

end # module SaturatedTaylor
#/ End module
