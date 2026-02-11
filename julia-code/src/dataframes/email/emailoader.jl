#= Module to load email dataset

=#
#/ Start module
module EmailLoader

#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase

#/ Modules, directories
import Meris.EMAILDIR as EMAILDIR

#################
### FUNCTIONS ###
function load(;
    DIR=GOWALLADIR * "raw-data/email-Enron.txt.gz"
)

    df = CSV.read(DIR, DataFrame, header=[:user_id, :time_z, :lat, :long, :location_id])

    # Consider dates in format YYY-MM-DD as samples and POIs as components.
    # Counts are the number of people entered in that POI in that day.
    
    dates = [row.time_z[1:10] for row in eachrow(df)]
    df.time = dates
    sdf = select(df, [:user_id, :location_id, :time])
    
    sdf = transform(
      groupby(sdf, [:time, :location_id]),
      :user_id => length => :counts
    )
    
    sdf = transform(
      groupby(sdf, [:time]),
      :counts => sum => :nreads
    )
    
    rename!(sdf, :location_id => :component_id, :time => :sample_id)
    sdf.class .= "gowalla"
    select!(sdf, [:class, :component_id, :sample_id, :counts, :nreads])

    return sdf
end

end # module GowallaLoader
#/ End module

