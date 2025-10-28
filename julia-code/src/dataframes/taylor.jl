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
function compute(df::DataFrame, idcolname; minoccupancy::Float64=1e-2, maxfrequency::Float64=1e-2)
    #~ Compute the (log) relative frequency of each of the "species"
    @transform!(df, :frequency = :counts ./ :nreads)
    @transform!(df, :logfrequency = log.(:frequency))
    #~ Compute some summary statistics
    nsamples = length(unique(df[!,:sample_id]))
    minoccupancy = minoccupancy < 1 ./ nsamples ? 1 ./ nsamples : minoccupancy
    sdf = @chain df begin
        @by(
            idcolname,
            :noccurences = length($idcolname),
            :occupancy = length($(idcolname)) ./ nsamples,
            :meanfrequency = mean(:frequency),
            :varfrequency = var(:frequency, corrected=false),
        )
        #~ Omit those with an occupancy below a specified threshold
        #~ note: if the threshold is zero, then omit hapax legomenas
        #        (those that occur only a single time within the entire document)
        @subset(:occupancy .> minoccupancy, :varfrequency .> 0.)
        #~ Select only those with a rel. frequency below the given maximum, and positive variance
        #  note: The first filters out those where, roughly, a Poisson distribution approximates
        #        a Binomial distribution well.
        @subset(:meanfrequency .< maxfrequency)
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
