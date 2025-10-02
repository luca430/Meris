#= Simple module to compute the AFD given a DataFrame =#
#/ Start module
module AFD

#/ Packages
using DataFrames, DataFramesMeta
using StatsBase

#################
### FUNCTIONS ###
function compute(df::DataFrame)
    #~ Filter
    
	  #~ Compute the (log) relative frequency of each of the "species"
    @transform!(df, :frequency = :counts ./ :nreads)
    @transform!(df, :logfrequency = log.(:frequency))
    #~ Compute some summary statistics
    nsamples = length(unique(df[!,:sample_id]))
    sdf = @chain df begin
        @by(
            :species_id,
            :occupancy = length(:sample_id) ./ nsamples,
            :meanlog = mean(:logfrequency),
            :stdlog  = std(:logfrequency, corrected=false)
        )
    end
    df = innerjoin(df, sdf, on=:species_id)
    df = @chain df begin
        @subset(:occupancy .≈ 1.0)
	      @transform(:z = (:logfrequency .- :meanlog) ./ :stdlog)
        @select(:sample_id,:species_id,:z)
    end
    return df
end

end # module AFD

#/ End module
