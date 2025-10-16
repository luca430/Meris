#= Simple module to compute the AFD given a DataFrame =#
#/ Start module
module AFD

#/ Packages
using DataFrames, DataFramesMeta
using StatsBase

#################
### FUNCTIONS ###
function compute(df::DataFrame, idcolname::Symbol; minoccupancy::Float64=1e-1)
	  #~ Compute the (log) relative frequency of each of the "species"/"component"
    @transform!(df, :frequency = :counts ./ :nreads)
    @transform!(df, :logfrequency = log.(:frequency))
    #~ Compute some summary statistics
    nsamples = length(unique(df[!,:sample_id]))
    sdf = @chain df begin
        @by(
            idcolname,
            :occupancy = length($(idcolname)) ./ nsamples,
            :meanlog = mean(:logfrequency),
            :varlog  = var(:logfrequency, corrected=false),
            # :stdlog  = std(:logfrequency, corrected=false)
        )
        #~ Take the occupation number into account
        #~ this means that μ → o⋅μ and σ² → o⋅[σ²+μ²(1-o)], where o the occupancy
        @transform(:meanlog = :meanlog .* :occupancy)
        @transform(:varlog = :varlog .+ :meanlog.^2 .* (1 .- :occupancy))
        @transform(:varlog = :varlog .* :occupancy)
        #~ Perform a log-transform on the mean-frequency (needed for lognormal)
        @subset(:varlog .> 0.0)
    end
    df = innerjoin(df, sdf, on=idcolname)
    df = @chain df begin
        @subset(:occupancy .> minoccupancy)
	      @transform(:z = (:logfrequency .- :meanlog) ./ sqrt.(:varlog))
        @select(:sample_id,$(idcolname),:z)
    end
    return df
end

end # module AFD

#/ End module
