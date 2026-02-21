#= Module to load stock volume data =#
#/ Start module
module FinanceLoader

#/ Packages
using Glob
using CSV, DataFrames, DataFramesMeta

#/ Modules, directories
import Meris.FINANCEDIR as FINANCEDIR

#################
### FUNCTIONS ###
"Load all financial stock volume data into a single DataFrame"
function load(
    ;
    DIR = FINANCEDIR * "raw-data/",
    filterdata    = true,
    minsamples    = 30,
    minreads      = 100_000,
    mincomponents = 100
    )
    #~ Gather list of files
    files = filter(f -> endswith(f, ".csv"), readdir(DIR, join=true))
    dfs = DataFrame[]

    #~ Read all markets into a single DataFrame
    for f in files
        tmp = CSV.read(f, DataFrame)
        # use filename (without extension) as class
        tmp.class .= split(splitext(basename(f))[1], "-")[1]
        push!(dfs, tmp)
    end
    df = vcat(dfs...)
    
    #~ Rename for consistency
    rename!(df, :ticker => :component_id, :total_volume => :counts)
    #~ Create DataFrame
    df = @chain df begin
        @subset(:counts .> 0)
        @groupby(:sample_id)
        @combine(:class, :sample_id, :component_id, :counts, :nreads = sum(:counts))
    end

    if filterdata
        #~ filter data
        @subset!(df, :nreads .> minreads)
        summarydf = @chain df begin
            @by(
                :class,
                :nsamples = length(:sample_id),
                :ncomponents = length(unique(:component_id))
            )
            @subset(:nsamples .> minsamples, :ncomponents .> mincomponents)
        end
        @subset!(df, :class .∈ Ref(summarydf.class))
    end
    #~ Return
    return df
end

end # module FinanceLoader
#/ End module

