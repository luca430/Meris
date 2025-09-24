#= Simple module to sample from books
   Books are 'special' in a way, as there is no noise, the complete size of the universe is
   known a priori, and samples of any size smaller than the universe can be taken.
=#
#/ Start module
module BookSampler

#/ Packages
using CSV
using DataFrames
using Random
using StatsBase

#/ Modules
import Moira.BOOKDIR as BOOKDIR

#################
### FUNCTIONS ###
"Take samples from a book and compute vocabulary size"
function samplevocabsize(
    Nv::Vector{Int};
    n::Int = 32,
    rng = Random.Xoshiro(42*n),
    DIR = BOOKDIR * "chinese/",
    FILENAME = "utf-story-of-the-stone.txt"
)
    words = load_book(; DIR=DIR, FILENAME=FILENAME)
    #~ For each N in Nv, sample N words n times and compute the vocabulary size V(N)
    V = zeros(Int, length(Nv), n)
    for i in eachindex(Nv), k in 1:n
        V[i,k] = _sample(words, Nv[i], rng=rng)
    end
    return V
end

########################
### HELPER FUNCTIONS ###
function load_book(;
    DIR = BOOKDIR * "chinese/",
    FILENAME = "utf-story-of-the-stone.txt"
)
    #/ Load into DataFrame
    df = CSV.read(DIR*FILENAME, DataFrame, header=["utfwords"])
    #~ Convert default String7 to String
    transform!(df, :utfwords => ByRow(String) => :utfwords)
    #~ Return a vector of String(s)
    return df.utfwords
end

"""
    _sample

Take a sample from a vector of 'words'. Words can be anything, as long as their String is unique.
"""
function _sample(
    words::Vector{String}, N::Int;
    rng = Random.Xoshiro(42)
)
    #~ Sample and count the vocabulary size V
    s = sample(rng, words, N, replace=false)
    return length(unique(s))
end
    
end # module BookSampler
#/ End module
