#= Taylor's sampler =#
#/ Start module
module HeapsGenerator

#/ Packages
using Distributions
using StatsBase
using Random
using SparseArrays
using JLD2

import Meris.MDistributions as MDist
import Meris.DATADIR as DATADIR
HEAPSDIR = DATADIR * "heaps/synthetic/"

#################
### FUNCTIONS ###
"""
    observe

Observe system by multinomial sampling out of S components, each with prob. `p[i] = θ[i]/Θ`, where
`Θ = sum(θ)`, and `θ ~ TemperedPareto(γ,1/φ)` with `φ` the characteristic scale of the system.
, when `N≫φ` 
"""
function generate(;
    N::Int=10^4,
    K::Int=10^3,
    γ::Float64=0.5,
    ε::Float64=1e0,
    S::Int=10^5,
    φ::Float64=1e2*S,
    nobservationslengths::Int = 31,
    rng=Random.Xoshiro(43*N)
)
    #~ Generate random propensities
    TPareto = MDist.TemperedPareto(γ, 1/φ, ε)
    θ = MDist.rand(rng, TPareto, S)
    #~ Specify noise distribution
    Pξ = Distributions.Gamma()

    #~ Log-divide N and count vocabulary size for each N'
    N = floor.(Int, exp10.(range(0, log10(N), nobservationslengths))) |> unique
    V = zeros(Int, length(N))
    counts = zeros(Int, S, K)
    for i in eachindex(N)
        ΔN = i > 1 ? N[i] - N[i-1] : 1
        for k in 1:K
            ξ = rand(rng, Pξ, S)
            p = θ .* ξ / sum(θ .* ξ)
            #~ Generate ΔN new counts
            counts[:,k] = counts[:,k] .+ rand.(rng, Poisson.(ΔN.*p))
            #~ Compute vocabulary size
            V[i] += count(x -> x > 0, counts[:,k])
        end
    end
    
    return (; N=N, V=V, params=(; N=N, K=K, S=S, φ=φ, γ=γ, ε=ε))
end

#############################
### DATA HELPER FUNCTIONS ###
function save_heaps(result; filename="synthetic-heaps.jld2", DIR=HEAPSDIR)
    mkpath(DIR)
    jldsave(DIR*filename, N=result.N, V=result.V, params=result.params)
end

end # module TaylorSampler
#/ End module
