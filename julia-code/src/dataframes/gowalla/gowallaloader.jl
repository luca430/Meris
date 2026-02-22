#= Module to load gowalla dataset

=#
#/ Start module
module GowallaLoader

#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase

using Meris

#/ Modules, directories
import Meris.GOWALLADIR as GOWALLADIR

#################
### FUNCTIONS ###
"Load all papers, put them into a single DataFrame"
function load(
    ;
    DIR = GOWALLADIR * "raw-data/",
    FILENAME = "loc-gowalla_totalCheckins.txt.gz",
    filterdata=true,
    minreads::Int=100_000,
    mincomponents::Int=1_000,
    minsamplecomponents::Int=500,
    minsamples::Int=30,         
)
    df = CSV.read(
        DIR*FILENAME, DataFrame, header=[:user_id, :time_z, :lat, :long, :location_id]
    )

    # Consider dates in format YYY-MM-DD as samples and POIs as components.
    # Counts are the number of people entered in that POI in that day.    
    dates = [row.time_z[1:10] for row in eachrow(df)]
    df.time = dates
    df = select(df, [:user_id, :location_id, :time])
    
    df = transform(
      groupby(df, [:time, :location_id]),
      :user_id => length => :counts
    )
    
    df = transform(
      groupby(df, [:time]),
      :counts => sum => :nreads
    )
    
    rename!(df, :location_id => :component_id, :time => :sample_id)
    df.class .= "gowalla"
    select!(df, [:class, :component_id, :sample_id, :counts, :nreads])

    if filterdata
        #~ filter data
        df = Meris.DataTools.df_filter(
            df;
            minreads=minreads,
            mincomponents=mincomponents,
            minsamplecomponents=minsamplecomponents,
            minsamples=minsamples
        )
    end

    return df
end

end # module GowallaLoader
#/ End module

