#= Simple module with functions to illustrate saturation of Heap's law when a specific
   characteristic scale of the system φ has been surpassed.
=#
#/ Start module
module SaturatedHeaps

#/ Packages
using StatsBase
using Random

#################
### FUNCTIONS ###
function compute_pseudocounts(φ::Float64, γ::Float64; S::Int = 2^15)
    α = map(i -> i^(-γ), 1:S)
    Z = sum(α)
    α = (φ / Z) .* α
    p = α ./ φ
    return (α, p)
end

function compute_latentvariables(α::Array{Float64})
	  φ = sum(α)
    return α ./ φ
end

# function compute_latentvariables(γ::Float64; S::Int = 256)
#     α = compute_pseudocounts(γ; S=S)
# 	  φ = sum(α)
#     return α ./ φ
# end

function generate_observation(N::Int, p::Array{Float64}; rng=Random.Xoshiro(42*N))
    S = length(p)
	  n = StatsBase.sample(rng, 1:S, N; replace=true)
    return n
end

function compute_vocabularysize(observation::Array)
	  cm = StatsBase.countmap(observation)
    return length(unique(keys(cm)))
end

########################
### HELPER FUNCTIONS ###
function compute_heaps(φ::Float64, γ::Float64; nobservations::Int = 50, s::Int = 23)
    #~ Pre-compute pseudo-counts and relative frequencies
    α, p = compute_pseudocounts(φ, γ)
    φ = sum(α)
    # p = compute_latentvariables(α)
    # idx = findlast(x -> x < 1, φ.*p)
    
	  #~ Distribution observation lengths N logarithmically around φ
    Nmin = 10
    Nmax = 100φ
    N = floor.(Int, exp10.(range(log10(Nmin), log10(Nmax), s)))
    V = zeros(nobservations, s)
    rng = Random.Xoshiro(nobservations*s)

    labels = collect(1:length(p))

    for k in 1:nobservations, i in 1:s
        # obs = generate_observation(N[i], p; rng=rng)
        obs = StatsBase.sample(rng, labels, Weights(p), N[i]; replace=true)
        V[k,i] = length(unique(obs))
    end
    V = dropdims(mean(V, dims=1), dims=1) ./ length(α)
    return N, V 
end

end # module SaturatedHeaps
#/ End module
