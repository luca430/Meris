#= Simple module to compute the mean and frequency for a given DataFrame
   Useful for investigating Taylor's law
=#
#/ Start module
module Taylor

#/ Packages
using DataFrames, DataFramesMeta
using StatsBase

#################
### FUNCTIONS ###
""
function compute(df::DataFrame, idcolname; minoccupancy::Float64=1e-1)
    #~ Compute the (log) relative frequency of each of the "species"
    @transform!(df, :frequency = :counts ./ :nreads)
    @transform!(df, :logfrequency = log.(:frequency))
    #~ Compute some summary statistics
    nsamples = length(unique(df[!,:sample_id]))
    sdf = @chain df begin
        @by(
            idcolname,
            :occupancy = length($(idcolname)) ./ nsamples,
            :meanfrequency = mean(:frequency),
            :varfrequency = var(:frequency, corrected=false),
            :meanlog = mean(:logfrequency),
            :varlog  = var(:logfrequency, corrected=false)
        )
        
        #~ Take the occupation number into account
        #~ this means that μ → o⋅μ and σ² → o⋅[σ²+μ²(1-o)], where o the occupancy
        @transform(:meanfrequency = :meanfrequency .* :occupancy)
        @transform(:varfrequency = :varfrequency .+ :meanfrequency.^2 .* (1 .- :occupancy))
        @transform(:varfrequency = :varfrequency .* :occupancy)
        # @transform(:meanlog = :meanlog .* :occupancy)
        # @transform(:varlog = :varlog .+ :meanlog.^2 .* (1 .- :occupancy))
        # @transform(:varlog = :varlog .* :occupancy)
        #~ Perform a log-transform on the mean-frequency (needed for lognormal)
        @subset(:occupancy .> minoccupancy, :varfrequency .> 0.0)
    end
    return sdf
end

end # module Taylor
#/ End module
