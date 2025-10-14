#= Module to sample from the LEGO dataset
   LEGO set can be obtained from https://rebrickable.com/downloads/
   Of particular interest to our use-case is the `inventory_parts`.
=#
#/ Start module
module LegoSampler

#/ Packages
using CSV
using DataFrames, DataFramesMeta
using Random
using StatsBase

#/ Modules, directories
import Meris.LEGODIR as LEGODIR

#################
### FUNCTIONS ###
"Take n samples of sizes Nv=[N1,N2,...] from the full LEGO catalogus and compute vocabulary size"
function samplevocabsize(
    Nv::Vector{Int};
    n::Int = 32,
    bagoflegos::Union{DataFrame,Nothing}=nothing,
    minquantity::Int = 100,      #~ Min. no of LEGO pieces in a set
    mindistinctpieces::Int = 30, #~ Min. no of *distinct* LEGO pieces in a set
    rng = Random.Xoshiro(42*minquantity*mindistinctpieces),
    DIR = LEGODIR,
    FILENAME = "inventory_parts.csv"
)
    if isnothing(bagoflegos)
        ldf = filterlegos(;
            minquantity=minquantity, mindistinctpieces=mindistinctpieces,
            DIR=DIR, FILENAME=FILENAME
        )
        bagoflegos = @select(ldf, :species_id, :counts)
    end
    #~ For each N in Nv, sample N words n times and compute the vocabulary size V(N)
    V = zeros(Int, length(Nv), n)
    for i in eachindex(Nv), k in 1:n
        V[i,k] = _samplevocabsize(Nv[i]; bagoflegos=bagoflegos, rng=rng)
    end
    return V
end


"""
    computevocabsize

Compute vocabulary size of LEGO sets of sufficient size
"""
function computevocabsize(;
    minquantity::Int = 50,       #~ Min. no of LEGO pieces in a set
    mindistinctpieces::Int = 50, #~ Min. no of *distinct* LEGO pieces in a set
    aggregate = false,
    returnsummary = true,
    DIR = LEGODIR,
    FILENAME = "inventory_parts.csv"
)
	  #/ Load into DataFrame
    #~ When `returnsummary=true`, return the summarizing statistics
    legodf = filterlegos(;
        minquantity=minquantity, mindistinctpieces=mindistinctpieces,
        DIR=DIR, FILENAME=FILENAME, returnsummary=returnsummary
    )
    #/ Aggregate the data by computing, for each unique sample size, the mean vocabularysize.
    #  This is useful for investigating Heap's law.
    if aggregate
        sdf = @by(legodf, :documentsize, :meanvocabularysize = mean(:vocabularysize))
        return sdf
    end
    return legodf
end

########################
### HELPER FUNCTIONS ###
"""
    filterlegos

Filter the LEGO dataset by including only sets with sufficient subvocabulary size and total
number of blocks
"""
function filterlegos(;
    ldf::Union{DataFrame,Nothing}=nothing,
    minquantity = 100,
    mindistinctpieces = 50,
    renamecols = true,
    returnsummary = false,
    DIR = LEGODIR,
    FILENAME = "inventory_parts.csv"
)
	  if isnothing(ldf)
        ldf = CSV.read(DIR*FILENAME, DataFrame)
    end
    #/ Select only relevant columns
    ldf = @select(ldf, :inventory_id, :part_num, :quantity)
    #~ Convert default String7 or String31 to String
    transform!(ldf, :part_num => ByRow(String) => :part_num)
    #/ Compute summary statistics and use it to filter
    sdf = @chain ldf begin
        @groupby(:inventory_id)
        @combine(:totalquantity = sum(:quantity), :distinctpieces = length(unique(:part_num)))
        @subset(:totalquantity .> minquantity, :distinctpieces .> mindistinctpieces)
    end
    if returnsummary
        #~ Rename columns for consistency
        sdf = @chain sdf begin
	          @rename(
                :documentsize = :totalquantity,
                :vocabularysize = :distinctpieces
            )
            @select(:documentsize, :vocabularysize)
        end
        return sdf
    end
    
    fdf = innerjoin(ldf, sdf, on=:inventory_id)
    
    if renamecols
        #~ Rename columns for consistency
        fdf = @chain fdf begin
            @rename(
	              :sample_id = :inventory_id,
                :species_id = :part_num,
                :counts = :quantity,
                :reads = :totalquantity,
                :vocabularysize = :distinctpieces
            )
            #~ [for now, omit vocabularysize size as it is not needed]
            @select(Not(:vocabularysize))
        end
    end
    return fdf
end

"""
    _samplevocabsize

Take sample of size `N` from the full LEGO catalogus and compute vocabulary size
"""
function _samplevocabsize(
    N::Int;
    bagoflegos::Union{DataFrame,Nothing}=nothing,
    minquantity::Int = 100,      #~ Min. no of LEGO pieces in a set
    mindistinctpieces::Int = 30, #~ Min. no of *distinct* LEGO pieces in a set
    rng = Random.Xoshiro(42*minquantity*mindistinctpieces),
    DIR = LEGODIR,
    FILENAME = "inventory_parts.csv"
)
    if isnothing(bagoflegos)
        ldf = filterlegos(;
            minquantity=minquantity, mindistinctpieces=mindistinctpieces,
            DIR=DIR, FILENAME=FILENAME
        )
        bagoflegos = @select(ldf, :species_id, :counts)
    end
    s = _sample(bagoflegos, N, rng=rng)
    V = length(unique(s))
    return V
end

    
"""
    _sample

Take a sample from a DataFrame of LEGO pieces. Use their counts as weights.
"""
function _sample(
    legos::DataFrame, N::Int;
    rng = Random.Xoshiro(42)
)    
    #~ Compute weights, and return sample from catalogus
    w = Weights(legos[!,:counts])
    s = sample(rng, legos[:,:species_id], w, N, replace=false)
    return s
end

##################################
### DATA ACQUISITION FUNCTIONS ###
"Download relevant LEGO dataset from https://rebrickable.com/downloads/"
function download(;
    URL = "https://cdn.rebrickable.com/media/downloads/inventory_parts.csv.zip?1758697954.19653",
    OUTDIR = LEGODIR
)
    mkpath(OUTDIR)
	  zpath = joinpath(OUTDIR, "inventory_parts.zip")
    run(`wget -O $(zpath) $(URL)`)
    run(`unzip -o $(zpath) -d $(OUTDIR)`)
    nothing
end



end # module LegoSampler
#/ End module
