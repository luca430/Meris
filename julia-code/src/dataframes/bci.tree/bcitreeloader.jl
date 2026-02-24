#= Module to load and sample from the Barro Colorade dataset
   Raw data can be obtained from https://doi.org/10.15146/5xcp-0d46, which lists `.RData` files
   that contain the information that is needed.
   Relevant `.RData` files are
   > `bci.tree.zip`      [contains `bc.tree[1-8].RData` files which contains indiv. tree records]
   > `bci.sptable.rdata` [contains meta data to unique identify trees, if desired]
=#
#/ Start module
module BCITreeLoader

#/ Packages
using CSV
using RData
using DataFrames, DataFramesMeta
using Random
using StatsBase

using Meris

#/ Modules, directories
import Meris.TREEDIR as TREEDIR

#################
### FUNCTIONS ###
#~ IMPORTANT: This function only works with steps=2, it needs to be debugged
function load(
        ;
        DIR=TREEDIR * "raw-data/",
        joinquadrats=true,
        steps::Int=30,
        filterdata=true,
        minreads::Int=5000,
        mincomponents::Int=200,
        minsamplecomponents::Int=100,
        minsamples::Int=30
    );
    treedf = DataFrame()
    census_ids = String[]
    filenames = readdir(DIR)
    for FILE in filenames
        #~ Check if numeric number in filename, otherwise skip
        (!any(isnumeric, FILE)) && (endswith(".rdata", FILE)) && (continue)
        df = load_rdata(; rdatafilename=FILE, joinquadrats=joinquadrats, steps=steps)
        #~ Add census_id [note: kept as a String]
        census_id = filter(isnumeric, FILE)
        push!(census_ids, census_id)
        df[!, :census_id] = fill(census_id, nrow(df))
        #~ Filter out those that have an ExactDate
        filter!(:ExactDate => d -> !ismissing(d), df)
        #~ Add the, to us, relevant columns to the pooled DataFrame        
        append!(
            treedf,
            @select(df, :census_id, :treeID, :sp, :status, :quadrat, :ExactDate, :date)
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
    for quadrat in unique(treedf[!, :quadrat])
        _df = filter(:quadrat => q -> q == quadrat, treedf)
        nreads = nrow(_df)
        cm = countmap(_df[!, :sp])
        #~ Add each entry of the census to the DataFrame
        for (tree, count) in pairs(cm)
            push!(countdf, [quadrat, tree, count, nreads], promote=true)
        end
    end

    (joinquadrats) && (countdf = join_quadrats(countdf; steps=steps))

    countdf.class .= "BCI"
    #~Filter entries if desired (since this dataset is small the filter can be done as a last step)
    if filterdata
        countdf = Meris.DataTools.df_filter(
            countdf,
            minreads=minreads,
            mincomponents=mincomponents,
            minsamplecomponents=minsamplecomponents,
            minsamples=minsamples
        )
    end
    return countdf
end

########################
### HELPER FUNCTIONS ###
"""
    load_rdata

Load RData file into a DataFrame
"""
function load_rdata(; rdatafilename="bci.tree1.rdata", joinquadrats=false, steps=2)
    df = RData.load(TREEDIR * "raw-data/" * rdatafilename)[splitext(rdatafilename)[begin]]
    df = @transform(df, :treeID = Int.(:treeID))
    # (joinquadrats) && (df = join_quadrats(df; steps=steps))
    return df
end

"""
    join_quadrats

Coarse grain by joining quadrats
"""
function join_quadrats(df::DataFrame; steps=2, rng=Random.Xoshiro(42 * steps))
    samples = unique(df.sample_id)
    samples = samples[randperm(rng, length(samples))]

    for n in 1:(length(samples) ÷ steps)
        _samples = samples[(n-1)*steps+1 : n*steps]

        mask = df.sample_id .∈ Ref(_samples)
        df.sample_id[mask] .= "Q$(n)"
    end

    df = @chain df begin
        @groupby(:sample_id, :component_id)
        @combine(:counts = sum(:counts), :nreads)
        @groupby(:sample_id)
        @combine(:component_id, :counts, :nreads = sum(:nreads))
    end

    return df
end

# function join_quadrats(df::DataFrame; steps=2, rng=Random.Xoshiro(42 * steps))
#     #/ Create summary
#     qdf = @chain df begin
#         #~ Count no. of trees in quadrat
#         @by(:quadrat, :ntrees = length(:sp))
#         @orderby(:ntrees)
#     end
#     #/ Remove quadrats with least no. of trees until a power of 2 is reached
#     s = floor(Int, log(nrow(qdf)) / log(2.))
#     qdf = qdf[nrow(qdf)-2^s:end, :]
#     #/ For each step, join two quadrats into a single one
#     quadratpairs = qdf[!, :quadrat]
#     for n in 1:steps
#         quadrats = randperm(rng, length(quadratpairs))
#         offset = length(quadrats) ÷ 2
#         quadratpairs = [[quadratpairs[quadrats[i]], quadratpairs[quadrats[i+offset]]] for i in 1:offset]
#     end
#     #~ Flatten the pairs (of pairs (...)) of quadrats
#     quadratgroups = map(q -> vcat(q...), quadratpairs)
#     quadratlabels = Dict(q => i for (i, group) in enumerate(quadratgroups) for q in group)
#     #/ Relabel all quadrats in the original DataFrame
#     df[!, :quadrat] = map(quadrat -> get(quadratlabels, quadrat, missing), df[!, :quadrat])
#     filter!(:quadrat => label -> !ismissing(label), df)
#     return df
# end

end # module BCTreeLoader
#/ End module
