#= Module to load gowalla dataset

=#
#/ Start module
module GowallaLoader

#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase

#/ Modules, directories
import ..DataTools: filterdata
import Meris.GOWALLADIR as GOWALLADIR

#################
### FUNCTIONS ###
"Load all papers, put them into a single DataFrame"
function load(
    ;
    DIR = GOWALLADIR * "raw-data/",
    FILENAME = "loc-gowalla_totalCheckins.txt.gz",
    minsamples    = 30,
    minreads      = 10_000,
    mincomponents = 100,
    applyfilter   = true,
    reorder       = true,
    top           = nothing,
)
    df = CSV.read(
        DIR*FILENAME, DataFrame, header=[:user_id, :time_z, :lat, :long, :location_id]
    )
    # Consider dates in format YYY-MM-DD as samples and POIs as components.
    # Counts are the number of people entered in that POI in that day.    
    dates = [row.time_z[1:10] for row in eachrow(df)]
    df.time = dates
    sdf = @chain df begin
        @select(:user_id, :location_id, :time)
        @groupby(:time, :location_id)
        @combine(:counts = length(:user_id))
        @groupby(:time)
        @transform(:nreads = sum(:counts))
    end
    
    #~ Standardize and rename
    sdf.class .= "gowalla"
    sdf = @rename(sdf, :sample_id = :time, :component_id = :location_id)
    @select!(sdf, :class, :sample_id, :component_id, :counts, :nreads)

    if applyfilter
        #~ filter data
        sdf = filterdata(
            sdf; minsamples=minsamples, minreads=minreads, mincomponents=mincomponents,
            reorder=reorder, top=top
        )
    end

    return sdf
end

end # module GowallaLoader
#/ End module

