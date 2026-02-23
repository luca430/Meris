#= Module to load stock volume data =#
#/ Start module
module FinanceLoader

#/ Packages
using Glob
using CSV, DataFrames, DataFramesMeta

using Meris

#/ Modules, directories
import Meris.FINANCEDIR as FINANCEDIR

#################
### FUNCTIONS ###
"Load all financial data, put them into a single DataFrame"
function load(
        ;
        DIR=FINANCEDIR * "raw-data/",
        filterdata=true,
        minreads::Int=10_000,
        mincomponents::Int=500,
        minsamplecomponents::Int=200,
        minsamples::Int=30,
    )
    files = filter(f -> endswith(f, ".csv"), readdir(DIR, join=true))
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
        @groupby(:sample_id)
        @combine(:class, :sample_id, :component_id, :counts, :nreads = sum(:counts))
    end
    
    if filterdata
        df = Meris.DataTools.df_filter(df,
            minreads=minreads,
            mincomponents=mincomponents,
            minsamplecomponents=minsamplecomponents,
            minsamples=minsamples,
        )
    end
    #~ Return
    return df
end

end # module FinanceLoader
#/ End module

