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
    n::Int=32,
    bagoflegos::Union{DataFrame,Nothing}=nothing,
    minquantity::Int=100,      #~ Min. no of LEGO pieces in a set
    mindistinctpieces::Int=30, #~ Min. no of *distinct* LEGO pieces in a set
    rng=Random.Xoshiro(42 * minquantity * mindistinctpieces),
    DIR=LEGODIR,
    FILENAME="inventory_parts.csv"
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
        V[i, k] = _samplevocabsize(Nv[i]; bagoflegos=bagoflegos, rng=rng)
    end
    return V
end


"""
    computevocabsize

Compute vocabulary size of LEGO sets of sufficient size
"""
function computevocabsize(;
    minquantity::Int=50,       #~ Min. no of LEGO pieces in a set
    mindistinctpieces::Int=50, #~ Min. no of *distinct* LEGO pieces in a set
    aggregate=false,
    returnsummary=true,
    DIR=LEGODIR,
    FILENAME="inventory_parts.csv"
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
    parse_themes

Load and filter LEGO data based on theme. Collects theme from metadata files and returns a
DataFrame where each LEGO piece additionally has a theme so that it can be filtered or selected.
"""
function parse_themes(;
    DIR=LEGODIR,
    SETFILE="sets.csv",
    INVENTORYSETFILE="inventories.csv",
    INVENTORYPARTSFILE="inventory_parts.csv",
    minquantity=64,
    mindistinctpieces=32,
    standardize=true,
    returnthemes=true
)
    #~ Load the DataFrames
    setdf = CSV.read(DIR * SETFILE, DataFrame)
    invdf = CSV.read(DIR * INVENTORYSETFILE, DataFrame)
    invpartdf = CSV.read(DIR * INVENTORYPARTSFILE, DataFrame)

    #~ Remove some columns that we don't need
    @select!(setdf, :set_num, :num_parts, :theme_id)
    @select(invdf, :id, :set_num)
    @select!(invpartdf, :inventory_id, :part_num, :quantity, :color_id)

    #/ All three DataFrames will be combined to add the appropriate metadata to each individual
    #  LEGO brick. The reason for this is that we can filter on specific themes in order to,
    #  potentially, "fix" the appropriate scale for a complex component system of a specific type.
    superdf = innerjoin(invpartdf, invdf, on=:inventory_id => :id)
    superdf = innerjoin(superdf, setdf, on=:set_num)
    #~ Omit entries with missing `set_num`, as those cannot be sorted or selected
    @subset!(superdf, map(x -> !ismissing(x), :set_num))

    #/ Omit and rename some columns when desired
    if standardize
        #~ Compute the total no. of bricks in each set in the inventory
        sdf = @chain superdf begin
            @groupby(:inventory_id)
            @combine(:nreads = sum(:quantity), :distinctpieces = length(unique(:part_num)))
            @subset(:nreads .> minquantity, :distinctpieces .> mindistinctpieces)
        end
        superdf = innerjoin(superdf, sdf, on=:inventory_id)
        #~ Rename some columns
        @rename!(superdf, :component_id = :part_num, :counts = :quantity, :sample_id = :inventory_id)
        #~ Pieces with the same component_id may have distinct colors, so make here a unique
        #  id that combines the component_id and the color_id
        @transform!(superdf, :component_id = :component_id .* "-" .* string.(:color_id))
        #~ Select only necessary columns
        @select!(superdf, :sample_id, :component_id, :counts, :nreads, :theme_id)
    end

    #/ Construct DataFrame with the no. of sets for each theme
    if returnthemes
        colname = standardize ? :sample_id : :inventory_id
        themedf = @by(superdf, :theme_id, :nsets = length(unique($(colname))))
        return superdf, themedf
    end
    return superdf
end

"""
    select_theme

Select specific theme.
If no theme is given (`theme_id=nothing`), selects the theme with the largest number of sets.
"""
function select_theme(df::DataFrame, themedf::DataFrame; theme_id=nothing)
    if isnothing(theme_id)
        mostcommontheme_id = themedf[!, :theme_id][argmax(themedf[!, :nsets])]
        return @subset(df, :theme_id .== mostcommontheme_id)
    end
    return @subset(df, :theme_id .== theme_id)
end

"""
    filterlegos

Filter the LEGO dataset by including only sets with sufficient subvocabulary size and total
number of blocks, omitting any information on type of the sets.
"""
function filterlegos(;
    ldf::Union{DataFrame,Nothing}=nothing,
    minquantity=100,
    mindistinctpieces=50,
    renamecols=true,
    returnsummary=false,
    DIR=LEGODIR,
    FILENAME="inventory_parts.csv"
)
    if isnothing(ldf)
        ldf = CSV.read(DIR * FILENAME, DataFrame)
    end
    #/ Select only relevant columns
    ldf = @select(ldf, :inventory_id, :part_num, :quantity, :color_id)
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
                :component_id = :part_num,
                :counts = :quantity,
                :nreads = :totalquantity,
                :vocabularysize = :distinctpieces
            )
            #~ Pieces with the same component_id may have distinct colors, so make here a unique
            #  id that combines the component_id and the color_id
            @transform(:component_id = :component_id .* "-" .* string.(:color_id))
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
    minquantity::Int=100,      #~ Min. no of LEGO pieces in a set
    mindistinctpieces::Int=30, #~ Min. no of *distinct* LEGO pieces in a set
    rng=Random.Xoshiro(42 * minquantity * mindistinctpieces),
    DIR=LEGODIR,
    FILENAME="inventory_parts.csv"
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

Take a random sample from a DataFrame of LEGO pieces. Use their counts as weights.
"""
function _sample(
    legos::DataFrame, N::Int;
    rng=Random.Xoshiro(42)
)
    #~ Compute weights, and return sample from catalogus
    w = Weights(legos[!, :counts])
    s = sample(rng, legos[:, :species_id], w, N, replace=false)
    return s
end

##################################
### DATA ACQUISITION FUNCTIONS ###
"Download relevant LEGO dataset from https://rebrickable.com/downloads/"
# [DEPRECATED]
# @TODO Ensure that the relevant files are of the same (or a specific) version
# For the relevant files, see `parse_themes`
function download(;
    URL="https://cdn.rebrickable.com/media/downloads/inventory_parts.csv.zip?1758697954.19653",
    OUTDIR=LEGODIR
)
    mkpath(OUTDIR)
    zpath = joinpath(OUTDIR, "inventory_parts.zip")
    run(`wget -O $(zpath) $(URL)`)
    run(`unzip -o $(zpath) -d $(OUTDIR)`)
    nothing
end



end # module LegoSampler
#/ End module
