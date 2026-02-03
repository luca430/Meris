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
    DIR=FINANCEDIR * "raw-data/",
    stopwords=true
)

    files = filter(f -> endswith(f, ".csv"), readdir(DIR, join=true))
    dfs = [CSV.read(f, DataFrame) for f in files]
    df = vcat(dfs...)

    return df
end

end # module FinanceLoader
#/ End module

