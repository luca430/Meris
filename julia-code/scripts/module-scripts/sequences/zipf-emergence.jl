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
function generate(
    N::Int;
    K::Int=512,
    α::Float64=1.0,
    γ::Float64=0.5,
    φ::Float64=1e5,
    ε::Float64=1e0,
    S::Int=2^12,
    rng=Random.Xoshiro(42*N)
)
    #~ Generate global component-specific intensities
    Pa = MDist.TemperedPareto(α, 1/φ, ε)
    a = MDist.rand(rng, Pa, S)
    
    c = 1e3 / S
    #~ Generate random propensities
    Pθ = MDist.TemperedPareto(γ, 1/φ/1e3, ε)
    #~ Specify noise distribution
    Pξ = Distributions.Gamma()
    
    #~ Allocate
    x = zeros(Int, S)

    for _ in 1:K
        #~ Sample probabilities
        ξ = rand(rng, Pξ, S)    #~ light-tailed perturbation [can be omitted]        
        θ = MDist.rand(rng, Pθ, S)        
        p = a .* θ .* ξ / sum(a .* θ .* ξ)
        #~ Multinomial sampling
        n = StatsBase.sample(rng, 1:S, Weights(p), N; replace=true)
        cmap = StatsBase.countmap(n)
        #~ Count
        for (label, counts) in cmap
            x[label] += counts
        end
    end

    #~ Compute Zipf [(rank, frequency)]
    Z = zipf(x)

    #~ Pre-compute points to make scatter plots faster
    k = length(Z.n)
    maxlength = 511
    step = floor(Int, k / maxlength)
    rp = vcat(collect(1:10), collect(11:k)[begin:step:end])
    np = Z.n[rp]
    
    return (; n=Z.n, rp=rp, np=np)
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
