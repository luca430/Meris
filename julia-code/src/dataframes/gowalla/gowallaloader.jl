#= Module to load gowalla dataset

=#
#/ Start module
module GowallaLoader

#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase

#/ Modules, directories
import Meris.GOWALLADIR as GOWALLADIR

#################
### FUNCTIONS ###
"Load all papers, put them into a single DataFrame"
function load(
    ;
    DIR = GOWALLADIR * "raw-data/",
    FILENAME = "loc-gowalla_totalCheckins.txt.gz",
    filterdata    = true,
    minsamples    = 30,
    minreads      = 10_000,
    mincomponents = 100          
)
    df = CSV.read(
        DIR*FILENAME, DataFrame, header=[:user_id, :time_z, :lat, :long, :location_id]
    )

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

    if filterdata
        #~ filter data
        @subset!(sdf, :nreads .> minreads)
        summarydf = @chain sdf begin
            @by(
                :class,
                :nsamples = length(:sample_id),
                :ncomponents = length(unique(:component_id))
            )
            @subset(:nsamples .> minsamples, :ncomponents .> mincomponents)
        end
        @subset!(sdf, :class .∈ Ref(summarydf.class))
    end

    return sdf
end

end # module GowallaLoader
#/ End module

