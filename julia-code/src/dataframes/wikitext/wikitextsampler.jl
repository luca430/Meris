#= Module to sample from the Wikitext-2 dataset
   Raw dataset can be obtained from;
   https://huggingface.co/datasets/Salesforce/wikitext/tree/main/wikitext-2-raw-v1
   which will be in a `.parquet` format. Parsing to raw word-per-text file is done in
   `parse-wikitext.py`, using the handy-dandy `pd.read_parquet`, as reading Parquet files
   in `Julia` is not super well-supported.
=#
#/ Start module
module WikitextSampler

#/ Packages
using CSV
using DataFrames, DataFramesMeta
using Random
using StatsBase

#/ Modules, directories
import Meris.WIKIDIR as WIKIDIR

#################
### FUNCTIONS ###
function computecdf(bagofwords::Vector{String}; ε::Float64 = 1e-12)
    cm = StatsBase.countmap(bagofwords)
    words = collect(keys(cm))
    n = length(words)
    frequencies = collect(values(cm))
    _order = sortperm(frequencies)
    permute!(words, _order)
    permute!(frequencies, _order)

    #~ Create list of unique frequencies
    #  but truncate it to machine precision
    uniquefrequencies = unique(frequencies)
    x = uniquefrequencies ./ n
    loweridx = searchsortedfirst(x, ε)
    uniquefrequencies = uniquefrequencies[loweridx:end]

    m = length(uniquefrequencies)
    F = Array{Float64}(undef, m)

    k = 1
    for i in eachindex(uniquefrequencies)
        while k <= n && frequencies[k] <= uniquefrequencies[i]
            k += 1
        end
        F[i] = (k - 1) / n
    end
    return (; F=F, t=uniquefrequencies)
end

function computecdf(;
    FILENAME = "wikitext-2-raw.txt",
    DIR = WIKIDIR,
)
    bagofwords = load_bagofwords(; FILENAME=FILENAME, DIR=DIR)
    return computecdf(bagofwords)
end

########################
### HELPER FUNCTIONS ###
function load_bagofwords(;
    FILENAME = "wikitext-2-raw.txt",
    DIR = WIKIDIR
)
    file = open(DIR*FILENAME)
    @info "name" DIR*FILENAME
    bagofwords = readlines(file)
    close(file)
    return bagofwords
end

end # module WikitextSampler
#/ End module
