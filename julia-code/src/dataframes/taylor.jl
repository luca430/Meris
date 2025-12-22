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
function compute(df::DataFrame, idcolname; minoccupancy::Float64=1e-2, maxfrequency::Float64=1e-2, occ::Bool=true)
        #~ Compute the (log) relative frequency of each of the "species"
        @transform!(df, :frequency = :counts ./ :nreads)
        # @transform!(df, :logfrequency = log.(:frequency))
        #~ Compute some summary statistics
        nsamples = length(unique(df[!,:sample_id]))
        minoccupancy = minoccupancy < 1 ./ nsamples ? 1 ./ nsamples : minoccupancy
        sdf = @chain df begin
                @by(
                        idcolname,
                        :noccurrences = length($idcolname),
                        :occupancy = length($(idcolname)) ./ nsamples,
                        :m = mean(:frequency),     #~ first moment
                        :m2 = mean(:frequency .^ 2),  #~ second moment
                        :m3 = mean(:frequency .^ 3),  #~ third moment
                        :m4 = mean(:frequency .^ 4),  #~ fourth moment
                        :s = var(:frequency)       #~ sample variance [unbiased]
                )
                @transform(:errorcov = (:noccurrences .- 1) .* :m3 ./ :noccurrences .^ 2)
                @transform(:varm = :s ./ :noccurrences)
                @transform(:vars = (:noccurrences .- 3) .* :s .^ 2)
                @transform(:vars = :vars ./ (:noccurrences .- 1))
                @transform(:vars = (:m4 .- :vars) ./ :noccurrences)
                @transform(:errorcorr = :errorcov ./ sqrt.(:varm .* :vars))
                #~ Omit those with an occupancy below a specified threshold
                #~ note: if the threshold is zero, then omit hapax legomenas
                #        (those that occur only a single time within the entire document)
                @subset(:occupancy .> minoccupancy, :s .> 0.)
                #~ Select only those with a rel. frequency below the given maximum, and positive variance
                #  note: The first filters out those where, roughly, a Poisson distribution approximates
                #        a Binomial distribution well.
                @subset(:m .< maxfrequency)
                #~ Take the occupation number into account
                #~ this means that μ → o⋅μ and σ² → o⋅[σ²+μ²(1-o)], where o the occupancy
                @transform(:varfrequency  = occ ? :occupancy .* (:s .+ :m.^2 .* (1 .- :occupancy)) : :s)
                @transform(:meanfrequency = occ ? :m .* :occupancy : :m)
                # @transform(:varfrequency = :s .* :occupancy)
                # @transform(:meanlog = :meanlog .* :occupancy)
                # @transform(:varlog = :varlog .+ :meanlog.^2 .* (1 .- :occupancy))
                # @transform(:varlog = :varlog .* :occupancy)
                #~ Perform a log-transform on the mean-frequency (needed for lognormal)
                # @rename(:meanfrequency = :m, :varfrequency = :s)
                @select(
                        :component_id,
                        :meanfrequency, :varfrequency,
                        :noccurrences, :occupancy,
                        :varm, :vars,
                        :errorcov, :errorcorr
                )
        end
        return sdf
end

end # module Taylor
#/ End module
