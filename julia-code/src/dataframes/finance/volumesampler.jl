#= Module to sample from Yahoo Finance volume data
=#
#/ Start module
module YahooSampler

#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase

import Meris.YAHOODIR as YAHOODIR

########################
### HELPER FUNCTIONS ###
function load_data(; filename = "yahoo-volumes.csv")
    return CSV.read(YAHOODIR * filename, DataFrame)    
end

end # module YahooSampler
#/ End module
