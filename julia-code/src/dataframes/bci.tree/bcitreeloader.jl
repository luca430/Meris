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
function load(
        ;
        DIR=TREEDIR * "raw-data/",
        joinquadrats=true,
        steps::Int=20,
        filterdata=true,
        minreads::Int=5000,
        mincomponents::Int=200,
        minsamples::Int=30,
        top=50
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
        countdf = Meris.DataTools.filterdata(
            countdf,
            minreads=minreads,
            mincomponents=mincomponents,
            minsamples=minsamples,
            top=top
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
    return df
end

"""
    join_quadrats

Coarse grain by joining quadrats
"""
function join_quadrats(df::DataFrame; steps=2, rng=Random.Xoshiro(42 * steps))
    (steps < 1) && throw(ArgumentError("steps must be >= 1"))

    samples = unique(df.sample_id)
    samples = samples[randperm(rng, length(samples))]

    n = 1
    i = 1
    while i <= length(samples)
        j = min(i + steps - 1, length(samples))
        _samples = samples[i:j]

        mask = df.sample_id .∈ Ref(_samples)
        df.sample_id[mask] .= "Q$(n)"

        n += 1
        i = j + 1
    end

    df = @chain df begin
        @groupby(:sample_id, :component_id)
        @combine(:counts = sum(:counts))
        @groupby(:sample_id)
        @combine(:component_id, :counts, :nreads = sum(:counts))
    end

    return df
end

end # module BCTreeLoader
#/ End module
