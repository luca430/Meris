#= Module to sample from the [static] Gutenberg dataset as
   provided by [...]
   https://zenodo.org/records/2422561
=#
#/ Start module
module GutenbergSampler

#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase

#/ Modules, directories
import Meris.GUTENBERGDIR as GUTENBERGDIR

#################
### FUNCTIONS ###

########################
### HELPER FUNCTIONS ###
function parse(; static=true)
end

end # module GutenbergSampler
#/ End module
