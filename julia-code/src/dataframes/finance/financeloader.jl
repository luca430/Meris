#= Module to load finance dataset

=#
#/ Start module
module FinanceLoader

#/ Packages
using Glob
using CSV, DataFrames, DataFramesMeta

#/ Modules, directories
import Meris.FINANCEDIR as FINANCEDIR

#################
### FUNCTIONS ###
"Load all financial data, put them into a single DataFrame"
function load(;
    DIR = FINANCEDIR * "raw-data/",
    stopwords = true
)

    files = filter(f -> endswith(f, ".csv"), readdir(DIR, join=true))

    dfs = DataFrame[]

    for f in files
        tmp = CSV.read(f, DataFrame)

        # use filename (without extension) as class
        tmp.class .= split(splitext(basename(f))[1], "-")[1]

        push!(dfs, tmp)
    end

    df = vcat(dfs...)

    rename!(df, :ticker => :component_id, :total_volume => :counts)

    df = @chain df begin
        @subset(:counts .> 0)
        @groupby(:sample_id)
        @combine(
            :class,
            :sample_id,
            :component_id,
            :counts,
            :nreads = sum(:counts)
        )
    end

    return df
end

end # module FinanceLoader
#/ End module

