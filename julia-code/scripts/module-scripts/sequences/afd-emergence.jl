#= Simple module to generate AFDs from latent variables =#
#/ Start module
module AFDGenerator

#/ Packages
using Distributions
using StatsBase
using Random
using SparseArrays
using FHist
using JLD2

import Meris.MDistributions as MDist
import Meris.DATADIR as DATADIR
AFDDIR = DATADIR * "afd/synthetic/"

#################
### FUNCTIONS ###
"""
    generate

Generate counts of components in a system by multinomial sampling out of S components, each with
prob. `p[i] = θ[i]/Θ`, where `Θ = sum(θ)`, and `θ ~ TemperedPareto(γ,1/φ)` with `φ` the
characteristic scale of the system.
"""
function generate(;
    N::Int=10^4,
    K::Int=10^4,
    γ::Float64=0.5,
    ε::Float64=1e0,
    S::Int=10^5,
    φ::Float64=1e1*S,
    Npcut::Float64=1000.0,
    rng=Random.Xoshiro(42*N)
)
    #~ Generate random propensities
    TPareto = MDist.TemperedPareto(γ, 1/φ, ε)
    θ = MDist.rand(rng, TPareto, S)
    #~ Specify noise distribution
    Pξ = Distributions.Gamma(1/(1-γ), 1/γ)

    #~ Allocate
    x = K > 2N ? spzeros(Int, K, S) : zeros(Int, K, S)
    pbar = zeros(Float64, S)
    
    for k in 1:K
        ξ = rand(rng, Pξ, S)
        p = θ .* ξ / sum(θ .* ξ)
        pbar .+= p / K
        x[k,:] .+= rand.(rng, Poisson.(N.*p))
    end
    
    #~ Compute occupancies of each component
    # o = count.(x -> x > 0, eachcol(x)) ./ K
    # commonidxs = findall(o -> o > 0.9, o)
    #  or, instead, select only those with Np̄<1
    commonidxs = findall(n -> n > Npcut, N*pbar)
    freqs = map(x -> x > zero(x) ? log(x) : zero(x), x[:,commonidxs] ./ N)
    μ = mean(freqs, dims=1)
    σ = var(freqs, dims=1)
    z = (freqs .- μ) ./ σ
    
    return (; z=z, params=(; N=N, Npcut=Npcut, K=K, S=S, φ=φ, γ=γ, ε=ε))
end

#############################
### DATA HELPER FUNCTIONS ###
function save_afd(result; filename="synthetic-afd.jld2", DIR=AFDDIR)
    mkpath(DIR)
    jldsave(DIR*filename, z=result.z, params=result.params)
end

end # module TaylorSampler
#/ End module
