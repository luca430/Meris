#= Simple module to sample from bag-of-words =#
#/ Start module
module WordSampler

#/ Packages
using CSV
using DataFrames
using Distributions
using Random
using StatsBase

#/ Modules, directories
import Meris.CORPUSDIR as CORPUSDIR

#################
### FUNCTIONS ###
"Take samples from a corpus (bag-of-words) and put them in a DataFrame"
function samplebagofwords(
    Nv::Vector{Int};
    replace = false,
    rng = Random.Xoshiro(42),
    DIR = CORPUSDIR * "brown/",
    FILENAME = "brown-corpus.txt"
)
    bagofwords = load_bagofwords(; DIR=DIR, FILENAME=FILENAME)
    #~ For each N in Nv, sample N words at random from the bag-of-words
    df = DataFrame(
        sample_id  = String[],
        species_id = String[],
        counts     = Int[],
        nreads     = Int[]
    )
    for i in eachindex(Nv)
        words = sample(rng, bagofwords, Nv[i], replace=replace)
        #~ Record the total no. of words (in this case, just Nv[i]), and the no. of counts
        #  for each distinct (unique) word
        cm = countmap(words)
        for (s, c) in cm
            push!(df, (i, s, c, Nv[i]), promote=true)
        end
    end
    return df
end

"""
    samplefrequencyofwords

Sample frequency of words obtained by sampling using the (relative) frequency of words in the bag.
First computes the relative frequencies, which are interpreted as probabilities. Then samples
either using a Multinomial [`replace=true`], or a Multivariate Hypergeometric [`replace=false]`.
"""
function samplefrequencyofwords(
    Nv::Vector{Int};
    replace = true,
    rng = Random.Xoshiro(42),
    DIR = CORPUSDIR * "brown/",
    FILENAME = "brown-corpus.txt"
)
    pdf = load_probabilities(; DIR=DIR, FILENAME=FILENAME)    
    (replace) && (return _sample_multinomial(Nv, pdf; rng=rng))
    return _sample_mvhypergeometric(Nv, pdf; rng=rng)
end

########################
### HELPER FUNCTIONS ###
"""
    _sample_multinomial
"""
function _sample_multinomial(
    N::AbstractVector{Int},
    pdf::DataFrame;
    rng = Random.Xoshiro(42)
)
	  #~ For each N in Nv, sample N words at random from a Multinomial with prob. vector p
    df = DataFrame(
        sample_id  = String[],
        species_id = String[],
        counts     = Int[],
        nreads     = Int[]
    )
    for i in eachindex(N)
        mult = Multinomial(N[i], pdf[!,:p])
        counts = rand(rng, mult)
        for k in eachindex(counts)
            if counts[k] > 0
                push!(df, (i, pdf[!,:species_id][k], counts[k], N[i]), promote=true)
            end
        end
    end
    return df
end

"""
    _sample_mvhypergeometric

Takes samples from multivariate hypergeometric distribution.
"""
function _sample_mvhypergeometric(N::AbstractVector{Int}, pdf::DataFrame; rng=Random.Xoshiro(42))
    #~ For each N in Nv, sample a total of N balls at random using a hypergeometric where the no.
    #  of successes (i.e., a ball of color s) is equal to the total no. of counts of that
    #  specific color
    df = DataFrame(
        sample_id  = String[],
        species_id = String[],
        counts     = Int[],
        nreads     = Int[]
    )
    S = size(pdf)[begin]    #~ Total no. of colors
    for i in eachindex(N)
        #~ Generate a random permutation of the words
        perm = randperm(rng, S)
        #~ Sample without replacement from that particular permutation
        counts = _sample_mvhypergeometric(pdf[!,:k][perm], N[i])
        for k in eachindex(counts)
            if counts[k] > 0
                idx = perm[k]
                push!(df, (i, pdf[!,:species_id][idx], counts[k], N[i]), promote=true)
            end
        end
    end
    return df
end

"""
    sample_mvhypergeometric

Use the hypergeometric approach to effectively draw N samples without replacement from a
bag-of-words with S categories each appearing n[s] times

!note: While using `sample(bagofwords, n, replace=false)` is probably much easier, this approach
       is [at least on paper] more efficient, and more closely resembles the "true" 
"""
function _sample_mvhypergeometric(
    n::AbstractVector{Int}, N::Int; rng=Random.Xoshiro(42*N)
)
    k = length(n) #~ No. of distinct colors
    x = zero(n)   #~ Stores no. of balls of color s
    M = sum(n)    #~ Total. no of balls currently in the urn
    
    # for i in eachindex(n)
    i = 0
    while i < k - 1 && N > 0
        i += 1
        if n[i] < M
            #/ Draw balls of color s following hypergeometric in N trials with
            #~ no. of successes: n[i]
            #~ no. of failures:  M - n[i]
            hypergeom = Hypergeometric(n[i], M - n[i], N)
            x[i] = rand(rng, hypergeom)
            M -= n[i] #~ Remove all balls of color s
            N -= x[i] #~ Subtract the no. of drawn balls
        else
            #/ If no more balls left over, all remaining balls are of color i
            x[i] = N
            N = 0
        end
    end
    #/ When there are still balls left over, these must be from the category k
    (i == k - 1) && (x[k] = N)
    return x
end

#############################
### DATA HELPER FUNCTIONS ###
function load_bagofwords(;
    DIR = CORPUSDIR * "brown/",
    FILENAME = "brown-corpus.txt"
)
    #/ Load into array
    bagofwords = CSV.read(DIR*FILENAME, DataFrame, header=["words"])
    #~ Convert to String
    transform!(bagofwords, :words => ByRow(String) => :words)
    return bagofwords.words
end

function load_probabilities(;
    DIR = CORPUSDIR * "brown/",
    FILENAME = "brown-corpus.txt"
)
	  #/ Load bag-of-words as array
    bagofwords = load_bagofwords(; DIR=DIR, FILENAME=FILENAME)
    #/ Compute relative frequencies
    cm = countmap(bagofwords)
    p  = values(cm)
    p  = p ./ sum(p)
    df = DataFrame(
        species_id = collect(keys(cm)),
        k = collect(values(cm)),
        p = p
    )
    return df
end

end # module WordSampler
#/ End module
