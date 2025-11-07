#= Module to parse and sample from the spatial patent dataset
   paper: https://doi.org/10.1038/s44260-025-00054-y
   data:  https://github.com/svalver/fractal-patenting/tree/main
=#
#/ Start module
module PatentSampler

using CSV
using DataFrames, DataFramesMeta

#/ Modules, directories
import Meris.PATENTDIR as PATENTDIR

#################
### FUNCTIONS ###
"Compute the survival function of the patent data across geographical locations"
function get_survival(df::DataFrame)
    x = sort(df[!,:count])
    C = cumsum(x)
    S = 1 .- C ./ C[end]
    return (; logx=log.(x), logS = log.(S))
end

########################
### HELPER FUNCTIONS ###

function load_data(;
    DIR = PATENTDIR,
    FILENAME = "full_biotech.csv"
)
    #~ Some entries in the data are `missing`, so silence warnings here
	  patentdf = CSV.read(DIR*FILENAME, DataFrame; normalizenames=true, silencewarnings=true)
    #~ drop missing entries
    dropmissing!(patentdf)
    #~ Create unique sample_id for each patent, extracted from (lon,lat)-coordinates
    @transform!(patentdf, :sample_id = hash.(string.(:longitude.*:latitude)))
    #~ Rename and select for consistency
    @rename!(patentdf, :component_id = :patent_id)
    @select!(patentdf, :sample_id, :component_id, :date, :country, :cpc)

    #~ For each sample (a geographical location, a city), count the no. of patents
    countdf = @chain patentdf begin
	      @by(
            :sample_id,
            :count = length(:component_id)
        )
        @subset(:count .> 32)
    end
    return countdf
end

end # module PatentSampler
#/ End module
