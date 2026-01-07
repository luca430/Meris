#= Module that contains methods to generate data from processes =#
#/ Start module
module PitmanYor

#/ Packages
using Random
using StatsBase
using AliasTables

#################
### FUNCTIONS ###
"""
    generate_samples

Generate samples from Pitman-Yor process with strength parameter θ>-α and discount parameter 0≤α<1.
When α=0, the process becomes a Dirichlet process.
"""
function generate_samples(
    N::Int; θ::Float64 = 10.0, α::Float64 = 0.1, rng=Random.Xoshiro(42*N))
    #/ Allocate
    counts = Int[]       #~ Counts of a specific category
    weights = Float64[]  #~ Weights for each category
    V = 0                #~ Total vocabulary size
    at = nothing         #~ AliasTable, needs to be `nothing` as it will be created conditionally

    for n in 1:N
        #~ Compute probability of sampling a new label
        Z = θ + n - 1
        p = (θ + α*V) / Z
        #~ Sample new or old label
        u = rand(rng)
        if u < p
            #~ Add new category and update counts
            push!(counts, 1)
            push!(weights, counts[end] - α)
            V += 1
            #~ Creation of AliasTable is O(v), with v the vocabulary size.
            #  So for large α
            at = AliasTable(weights)
        else
            #~ Sample existing category
            k = rand(rng, at)
            #~ Update counts, weights, and AliasTable
            counts[k] += 1
            weights[k] = counts[k] - α
            #~ in-place update of AliasTable
            AliasTables.set_weights!(at, weights)
        end
    end
    return at
end

"""
    generate_pseudocounts

Generate pseudocounts α by a Pitman-Yor process until a set of S distinct labels are assigned.
"""
function generate_pseudocounts(
    S::Int; θ::Float64 = 10.0, α::Float64 = 0.1, rng=Random.Xoshiro(42*S))
    #/ Allocate
    counts = Int[]       #~ Counts of a specific category
    weights = Float64[]  #~ Weights for each category
    V = 0                #~ Total vocabulary size
    at = nothing         #~ AliasTable, needs to be `nothing` as it will be created conditionally
    n = 0

    while length(counts) < S
        n += 1
        #~ Compute probability of sampling a new label
        Z = θ + n - 1
        p = (θ + α*V) / Z
        #~ Sample new or old label
        u = rand(rng)
        if u < p
            #~ Add new category and update counts
            push!(counts, 1)
            push!(weights, counts[end] - α)
            V += 1
            #~ Creation of AliasTable is O(v), with v the vocabulary size.
            #  So for large α
            at = AliasTable(weights)
        else
            #~ Sample existing category
            k = rand(rng, at)
            #~ Update counts, weights, and AliasTable
            counts[k] += 1
            weights[k] = counts[k] - α
            #~ in-place update of AliasTable
            AliasTables.set_weights!(at, weights)
        end
    end
    return at, counts
end

"""
!!!![work in progress]
    generate_concentrated_samples

Generates concentrated samples from Pitman-Yor process with strength parameter θ>-α and discount
parameter 0≤α<1, and an initial set of weights [concentration parameters] a.
Note that when α=0, the process becomes a Dirichlet process.
"""
function generate_samples(
    N::Int,
    wcounts::Array{Int64};
    w::Float64 = 0.2,
    θ::Float64 = 10.0,
    α::Float64 = 0.1,
    rng=Random.Xoshiro(42*N),
    return_counts::Bool = true
)
    #~ Pre-compute weights
    bweights = wcounts ./ sum(wcounts)
    A = sum(bweights)
    B = A * (1-w) / w
    #/ Allocate
    counts = Int[]
    weights = Float64[]
    V = 0
    
    at = nothing #AliasTable(weights)

    for n in 1:N
        #~ Compute probability of sampling a new label
        p = (θ + α*V + B) / (θ + n - 1 + B)
        #~ Sample new or old label
        u = rand(rng)
        if u < p
            #~ Add new category and update counts
            push!(counts, 1)
            push!(weights, counts[end] - α)
            V += 1
            #~ Creation of AliasTable is O(v), with v the vocabulary size.
            #  So for large α
            at = AliasTable(weights)
        else
            #~ Sample existing category
            k = rand(rng, at)
            #~ Update counts, weights, and AliasTable
            counts[k] += 1
            weights[k] = counts[k] - α
            #~ in-place update of AliasTable
            AliasTables.set_weights!(at, weights)
        end
    end
    (!return_counts) && (return at)
    return (at, counts)
end


########################
### HELPER FUNCTIONS ###
"""
		countvocabsize

Compute (count) vocabulary size of samples from a Pitman-Yor process.
Wrapper around `generate_pitmanyor_samples` to execute `n` processes for some values of N.
Multi-threading is enabled, as for large `N` sampling may take some time.
"""
function countvocabsize(
    N::Vector{Int};
    n = 30,
    θ=1.0,
    α=0.1,
    rng=Random.Xoshiro(42*n)
)
    V = zeros(Int, length(N), n)
    for i in eachindex(N)
        Threads.@threads for k in 1:n
            at = generate_samples(N[i]; θ=θ, α=α, rng=rng)
            V[i,k] += length(at)
        end
    end
    return V
end

end # module Process
#/ End module
