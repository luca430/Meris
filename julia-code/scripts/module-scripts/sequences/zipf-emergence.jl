#= Zipf's generator =#
#/ Start module
module ZipfGenerator

#/ Packages
using Distributions
using StatsBase
using JLD2
using Random

import Meris.MDistributions as MDist
import Meris.DATADIR as DATADIR
ZIPFDIR = DATADIR * "zipf/synthetic/"

#################
### FUNCTIONS ###
function generate(;
    N::Int=10^4,
    K::Int=10^3,
    γ::Float64=0.5,
    ε::Float64=1e0,
    S::Int=10^5,
    φ::Float64=1e2*S,
    rng=Random.Xoshiro(42*N)
)
    #~ Generate random propensities
    Pθ = MDist.TemperedPareto(γ, 1/φ, ε)
    #~ Specify noise distribution
    Pξ = Distributions.Gamma()
    
    #~ Allocate
    x = zeros(Int, S)

    θ = MDist.rand(rng, Pθ, S)
    for k in 1:K
        #~ Sample probabilities
        ξ = rand(rng, Pξ, S)    #~ light-tailed perturbation [can be omitted]
        p = θ .* ξ / sum(θ .* ξ)
        #~ Poisson sampling
        x .+= rand.(rng, Poisson.(N.*p))
        #~ Multinomial sampling
        # n = StatsBase.sample(rng, 1:S, Weights(p), N; replace=true)
        # cmap = StatsBase.countmap(n)
        # #~ Count
        # for (label, counts) in cmap
        #     x[label] += counts
        # end
    end

    #~ Compute Zipf [(rank, frequency)]
    Z = zipf(x)
    Snz = length(Z.n)

    #~ Pre-compute points to make scatter plots faster
    k = length(Z.n)
    _ranks = floor.(Int, exp10.(range(0, log10(Snz), 101))) |> unique
    _frequencies = Z.n[_ranks]
    
    return (; ranks=_ranks, freqs=_frequencies, params=(; N=N, K=K, S=S, φ=φ, γ=γ, ε=ε))
end

########################
### HELPER FUNCTIONS ###
function zipf(x::Array{T}) where {T<:Real}
    #~ Compute sorted counts
    xnz = x[x .> 0]
    ns = xnz ./ sum(xnz)
    ranks = sortperm(ns, rev=true)
    return (; n=ns[ranks])
end

#############################
### DATA HELPER FUNCTIONS ###
function save_zipf(result; filename="synthetic-zipf.jld2", DIR=ZIPFDIR)
    mkpath(DIR)
    jldsave(
        DIR*filename,
        ranks=result.ranks, freqs=result.freqs,
        params=result.params
    )
end

end # module ZipfGenerator
#/ End module
