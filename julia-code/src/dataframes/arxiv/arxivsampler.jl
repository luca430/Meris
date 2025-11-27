#= Module to investigate macro(-ecological) laws in arXiv paper

=#
#/ Start module
module arXivSampler

#/ Packages
using Glob
using CSV, DataFrames, DataFramesMeta
using Random, StatsBase

#/ Modules, directories
import Meris.ARXIVDIR as ARXIVDIR

#################
### FUNCTIONS ###
"Compute vocabulary size vs. total document size by appending documents in random order"
function computevocabsize(df::DataFrame; rng=Random.Xoshiro(42))
    #/ Construct vocabulary and dictionary
    #  note: dictionary is a set for quick comparison
    heapdf = DataFrame(documentsize=Int[], vocabularysize=Int[])
    vocabulary = String[]
    dictionary = Set{String}()
    #/ In random order, compute the vocabularysize for increasing document sizes
    sample_ids = unique(df[!,:sample_id])
    _order = randperm(rng, length(sample_ids))
    for id in sample_ids[_order]
        idxs = findall(df[!,:sample_id] .== id)
        for word in df[!,:component_id][idxs]
	          if !(word in dictionary)
                push!(vocabulary, word)
                push!(dictionary, word)
            end
        end
        
        _documentsize = isempty(heapdf[!,:documentsize]) ? length(idxs) :
                        last(heapdf[!,:documentsize]) + length(idxs)
        push!(heapdf, [_documentsize, length(vocabulary)])
    end
    return heapdf
end

function summarize(df::DataFrame)
	  sdf = @chain df begin
        @by(
            :sample_id,
            :totalcounts = length(:component_id),
            :vocabularysize  = length(unique(:component_id))
        )
    end
    return sdf
end

########################
### HELPER FUNCTIONS ###
"Load all papers of a specific category, put them into a single DataFrame"
function load_papers(;
    CATEGORY = "q-bio.PE",
    DIR = ARXIVDIR * "raw-text/"    
)
    #~ Specify directory
    SUBDIR = join(split(CATEGORY, "."), "/")
    _DIR = DIR * SUBDIR * "/"

    df = DataFrame(sample_id=String[], component_id=String[])

    for filename in readdir(_DIR)        
        paperdf = CSV.read(_DIR * filename, DataFrame, header=["component_id"])
        paperid = split(filename, ".")[begin]
        paperdf.sample_id = fill(paperid, nrow(paperdf))
        append!(df, paperdf)
    end
    return df
end

"Load *all* available papers of *all* categories"
function load_all(
    DIR = ARXIVDIR * "lemmatized-text/"
)
    df = DataFrame(sample_id=String[], component_id=String[])

    for (root, dirs, files) in walkdir(DIR)
        filenames = joinpath.(root, files)
        for (filename, paperid) in zip(filenames, files)
            if endswith(filename, ".txt")
                paperdf = CSV.read(filename, DataFrame, header=["component_id"])
                paperdf.sample_id = fill(paperid, nrow(paperdf))
                append!(df, paperdf)
            end
        end        
    end
    return df	  
end

end # module arXivSampler
#/ End module

