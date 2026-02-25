#= Module to load stock volume data =#
#/ Start module
module FinanceLoader

#/ Packages
using Glob
using CategoricalArrays
using CSV, DataFrames, DataFramesMeta

using Meris

#/ Modules, directories
# using Meris.DataTools
import ..DataTools: filterdata
import Meris.FINANCEDIR as FINANCEDIR


#################
### FUNCTIONS ###
"Load all financial data, put them into a single DataFrame"
function load(
    ;
    DIR = FINANCEDIR * "raw-data/",
    applyfilter   = true,
    minsamples    = 30,
    minreads      = 100_000,
    mincomponents = 100,
    resolution    = "daily",
    reorder       = true,
    top           = nothing
    )
    #~ Gather list of files
    files = filter(f -> endswith(f, "$(resolution)-volumes.csv"), readdir(DIR, join=true))
    dfs = DataFrame[]

    #~ Read all markets into a single DataFrame
    for f in files
        tmp = CSV.read(f, DataFrame)
        # use filename (without extension) as class
        class = split(splitext(basename(f))[1], "-")
        tmp.class .= class[1] .* "-" .* class[2]

        push!(dfs, tmp)
    end
    df = vcat(dfs...)
    
    #~ Rename for consistency
    rename!(df, :ticker => :component_id, :total_volume => :counts)
    #~ Create DataFrame
    df = @chain df begin
        @subset(:counts .> 0)
        @groupby(:class, :sample_id)
        @combine(:class, :sample_id, :component_id, :counts, :nreads = sum(:counts))
    end

    if applyfilter
        #~ filter data
        df = filterdata(
            df; minsamples=minsamples, minreads=minreads, mincomponents=mincomponents,
            reorder=reorder, top=top
        )
    end
    #~ Return
    return df
end

end # module FinanceLoader
#/ End module

