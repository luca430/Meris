#= Module to sample from the LEGO dataset
   LEGO set can be obtained from https://rebrickable.com/downloads/
   Of particular interest to our use-case are `inventory_parts.csv`,
   `inventories.csv`, and `sets.csv`.
=#
#/ Start module
module LEGOLoader

#/ Packages
using CSV
using DataFrames, DataFramesMeta
using Random
using StatsBase

using Meris

#/ Modules, directories
import ..DataTools: filterdata
import Meris.LEGODIR as LEGODIR

#######################
### FUNCTIONS ###
function load(
    ;
    DIR=LEGODIR * "raw-data/",
    SETFILE="sets.csv",
    INVENTORYSETFILE="inventories.csv",
    INVENTORYPARTSFILE="inventory_parts.csv",
    nthemes       = 20,
    minsamples    = 30,
    minreads      = 1000,
    mincomponents = 100,
    applyfilter   = true,
    reorder       = true,
    top           = nothing,        
    )

    df, themes_df = parse_themes(;
        DIR=DIR,
        SETFILE=SETFILE,
        INVENTORYSETFILE=INVENTORYSETFILE,
        INVENTORYPARTSFILE=INVENTORYPARTSFILE
        )

    rename!(df, :theme_id => :class)
    df.class .= string.(df.class)

    if applyfilter
        #~ filter data
        df = filterdata(
            df; minsamples=minsamples, minreads=minreads, mincomponents=mincomponents,
            reorder=reorder, top=top
        )
    end

    return df
    # if aggregate
    #     #~ aggregate as if it was a single `sample_id`
    #     aggdf = @chain df begin
    #         @groupby(:class, :sample_id, :component_id)
    #         @combine(:counts = sum(:counts))
    #     end
    #     @transform!(aggdf, :nreads = sum(:counts))
    #     return aggdf
    # end
    # return df
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
    INVENTORYPARTSFILE="inventory_parts.csv"
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

    #~ Compute the total no. of bricks in each set in the inventory
    sdf = @chain superdf begin
        @groupby(:inventory_id)
        @combine(:nreads = sum(:quantity), :distinctpieces = length(unique(:part_num)))
    end
    superdf = innerjoin(superdf, sdf, on=:inventory_id)
    #~ Rename some columns
    @rename!(superdf, :component_id = :part_num, :counts = :quantity, :sample_id = :inventory_id)
    #~ Pieces with the same component_id may have distinct colors, so make here a unique
    #  id that combines the component_id and the color_id
    @transform!(superdf, :component_id = :component_id .* "-" .* string.(:color_id))
    #~ Select only necessary columns
    @select!(superdf, :sample_id, :component_id, :counts, :nreads, :theme_id)

    #/ Construct DataFrame with the no. of sets for each theme
    themedf = @by(superdf, :theme_id, :nsets = length(unique(:sample_id)))
    return superdf, themedf
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

end # module LegoLoader
#/ End module
