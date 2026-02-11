#/ Start module
module FlightSampler

using CSV, DataFrames, DataFramesMeta

import Meris.FLIGHTDIR as FLIGHTDIR

#################
### FUNCTIONS ###
function load(
    ;
    FILENAME = "T_ONTIME_MARKETING.csv"
    )
    df = CSV.read(FLIGHTDIR*FILENAME, DataFrame)
    dropmissing!(df)
    return df
end


function count(df = nothing)
    (isnothing(df)) && (df = load())
    @rename!(df, :sample_id = :FL_DATE)
    #~ Reshape to long version
    airportdf = stack(
        df,
        [:ORIGIN_AIRPORT_ID, :DEST_AIRPORT_ID],
        variable_name = :type,
        value_name = :airport
    )
    flightdf = @chain airportdf begin
        @groupby(:sample_id, :airport)
        @combine(:counts = length(:airport))
    end
    return flightdf
    
    @transform!(df, :component_id = string.(:ORIGIN_AIRPORT_ID) .* string.(:DEST_AIRPORT_ID))
    
	  daydf = @chain df begin
        @groupby(:sample_id)
        @combine(:nreads = length(:sample_id))
    end
    #~ Count the invididual species per year per `speciesKey`
    countdf = @chain df begin
        @groupby(:sample_id, :component_id)
        @combine(:counts = length(:component_id))
    end
    return countdf
    #~ Rename some columns
    countdf = leftjoin(countdf, daydf, on=:sample_id)
    # @rename!(countdf, :sample_id = :year, :component_id = :speciesKey)
    return countdf
end

end # module FlightSampler
#/ End module
