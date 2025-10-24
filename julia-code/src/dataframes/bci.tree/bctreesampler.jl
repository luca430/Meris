#= Module to load and sample from the Barro Colorade dataset
   Raw data can be obtained from https://doi.org/10.15146/5xcp-0d46, which lists `.RData` files
   that contain the information that is needed.
   Relevant `.RData` files are
   > `bci.tree.zip`      [contains `bc.tree[1-8].RData` files which contains indiv. tree records]
   > `bci.sptable.rdata` [contains meta data to unique identify trees, if desired]
=#
#/ Start module
module BCITreeSampler

#/ Packages
using CSV
using RData
using DataFrames, DataFramesMeta
using Random
using StatsBase

#/ Modules, directories
import Meris.TREEDIR as TREEDIR
const RTREEDIR = TREEDIR * "RData/"
const CSVTREEDIR = TREEDIR * "csv/"

#################
### FUNCTIONS ###
function load_treedata(; filename=nothing, mincount=1, joinquadrats=true)
    treedf = DataFrame()
    census_ids = String[]
    filenames = isnothing(filename) ? readdir(RTREEDIR) : [filename]
    for FILE in filenames
        #~ Check if numeric number in filename, otherwise skip
        (!any(isnumeric, FILE)) && (continue)
        df = load_rdata(; rdatafilename=FILE, joinquadrats=joinquadrats)
        #~ Add census_id [note: kept as a String]
        census_id = filter(isnumeric, FILE)
        push!(census_ids, census_id)
        df[!,:census_id] = fill(census_id, nrow(df))
        #~ Filter out those that have an ExactDate
        filter!(:ExactDate => d -> !ismissing(d), df)
        #~ Record the total. no of tree (`nreads`) in that census, as it's needed later        
        # nreads = nrow(df)
        # df[!,:nreads] = fill(nreads, nreads)
        #~ Add the, to us, relevant columns to the pooled DataFrame        
        append!(
            treedf,
            @select(df, :census_id,:treeID,:sp,:status,:quadrat,:ExactDate,:date)
        )
    end
    #~ Filter out all those that have died in the meantime
    alivedf = @chain treedf begin
        @by(
            :treeID,
            :alive = all(:status .== "A") && length(:status) == length(census_ids)
        )
        @subset(:alive .== true)
    end
    treedf = innerjoin(treedf, alivedf, on=:treeID)
    #~ Allocate
    countdf = DataFrame(sample_id=String[], component_id=String[], counts=Int[], nreads=Int[])
    #~ Per plot [quadrat], count the no. of trees of each species
    for quadrat in unique(treedf[!,:quadrat])
        _df = filter(:quadrat => q -> q == quadrat, treedf)
        nreads = nrow(_df)
        cm = filter(x -> last(x) >= mincount, countmap(_df[!,:sp]))
        #~ Add each entry of the census to the DataFrame
        for (tree, count) in pairs(cm)
            (count > mincount) && (push!(countdf, [quadrat, tree, count, nreads], promote=true))
        end
    end    
    return countdf
end

########################
### HELPER FUNCTIONS ###
"""
    load_rdata

Load RData file into a DataFrame
"""
function load_rdata(; rdatafilename = "bci.tree1.rdata", joinquadrats=false)
	  df = RData.load(RTREEDIR*rdatafilename)[splitext(rdatafilename)[begin]]
    df = @transform(df, :treeID = Int.(:treeID))
    (joinquadrats) && (df = join_quadrats(df))
    return df
end

"""
    load_data

Load DataFrame from specific document
"""
function load_data(documentname; csvfilename = "raw-bci.tree-$(documentname).csv")
	  return CSV.read(CSVTREEDIR*csvfilename, DataFrame, delim=", ")
end

"""
    join_quadrats

Coarse grain by joining quadrats
"""
function join_quadrats(df::DataFrame; steps=2, rng=Random.Xoshiro(42*steps))
    #/ Create summary
    qdf = @chain df begin
        #~ Count no. of trees in quadrat
        @by(:quadrat, :ntrees = length(:sp))
        @orderby(:ntrees)
    end
    #/ Remove quadrats with least no. of trees until a power of 2 is reached
    s = floor(Int, log(nrow(qdf)) / log(2.))
    qdf = qdf[nrow(qdf)-2^s:end, :]
    #/ For each step, join two quadrats into a single one
    quadratpairs = qdf[!,:quadrat]
    for n in 1:steps
        quadrats = randperm(rng, length(quadratpairs))
        offset = length(quadrats) ÷ 2
        quadratpairs = [[quadratpairs[i], quadratpairs[i+offset]] for i in 1:offset]
    end
    #~ Flatten the pairs (of pairs (...)) of quadrats
    quadratgroups = map(q -> vcat(q...), quadratpairs)
    quadratlabels = Dict(q => i for (i,group) in enumerate(quadratgroups) for q in group)
    #/ Relabel all quadrats in the original DataFrame
    df[!,:quadrat] = map(quadrat -> get(quadratlabels, quadrat, missing), df[!,:quadrat])
    filter!(:quadrat => label -> !ismissing(label), df)
    return df
end

end # module BCTreeSampler
#/ End module
