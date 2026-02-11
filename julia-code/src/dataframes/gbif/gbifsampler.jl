#= 

Queries:
GBIF.org (11 February 2026) GBIF Occurrence Download  https://doi.org/10.15468/dl.cax2up
=#


#/ Start module
module GBIFSampler

using CSV, DataFrames, DataFramesMeta

import Meris.GBIFDIR as GBIFDIR

#################
### FUNCTIONS ###
function load(
    ;
    FILENAME="plantae-northamerica.csv",
    DIR = GBIFDIR,
    select=["speciesKey", "year"]
    )
    df = CSV.read(DIR*FILENAME, DataFrame; select=["basisOfRecord", "publishingOrgKey", "institutionCode", "collectionCode"])
    # dropmissing!(df)
    return df
end

"Get the counts for each year in the GBIF dataset"
function get_counts(; df = nothing)
    (isnothing(df)) && (df = load())
    #~ Get the total number of reads in a year
    #  [here, "reads" refer to the total number of occurences]
    readdf = @chain df begin
        @groupby(:year)
        @combine(:nreads = length(:speciesKey))
    end
    #~ Count the invididual species per year per `speciesKey`
    countdf = @chain df begin
        @groupby(:year, :speciesKey)
        @combine(:counts = length(:speciesKey))
    end
    #~ Rename some columns
    countdf = leftjoin(countdf, readdf, on=:year)
    @rename!(countdf, :sample_id = :year, :component_id = :speciesKey)
    return countdf
end

########################
### HELPER FUNCTIONS ###
"To reduce the file size, save just the relevant counts here"
function save_counts(
    countdf::DataFrame;
    DIR=GBIFDIR,
    FILENAME="plantae-northamerica-counts.csv"
    )
	  CSV.write(DIR*FILENAME, countdf)
end

end # module GBIFSampler
#/ End module
