#= Simple module to compute the AFD given a DataFrame =#
#/ Start module
module AFD

#/ Packages
using DataFrames, DataFramesMeta
using StatsBase

#################
### FUNCTIONS ###
function compute(
    df::DataFrame,
    idcolname::Symbol;
    maxfrequency::Float64=1e-2,
    minoccupancy::Float64=1e-2,
    occ::Bool=true,
    rescale_by_occupancy::Union{Nothing,Bool}=nothing,
    normalize_by_nreads::Bool=true
)
    occupancy_rescale = isnothing(rescale_by_occupancy) ? occ : rescale_by_occupancy
	  #~ Compute the (log) relative frequency/abundance of each "species"/"component"
    if normalize_by_nreads
        @transform!(df, :frequency = :counts ./ :nreads)
    else
        @transform!(df, :frequency = :counts)
    end
    @transform!(df, :logfrequency = log.(:frequency))
    #~ Compute some summary statistics
    nsamples = length(unique(df[!,:sample_id]))
    minoccupancy = minoccupancy < 1 ./ nsamples ? 1 ./ nsamples : minoccupancy
    sdf = @chain df begin
        @by(
            idcolname,
            :noccurences = length($(idcolname)),
            :occupancy = length($(idcolname)) ./ nsamples,
            :meanfrequency = mean(:frequency),
            :varfrequency = var(:frequency, corrected=false),
            :meanlog = mean(:logfrequency),
            :varlog  = var(:logfrequency, corrected=false),
        )
        #~ Omit those with an occupancy below a specified threshold
        #~ note: if the threshold is zero, then omit hapax legomenas
        #        (those that occur only a single time within the entire document)
        @subset(:occupancy .> minoccupancy, :varlog .> 0.)
        #~ Select only those with a rel. frequency below the given maximum, and positive variance
        #  note: The first filters out those where, roughly, a Poisson distribution approximates
        #        a Binomial distribution well.
        @subset(:meanfrequency .< maxfrequency)
        #~ Take the occupation number into account
        #~ this means that μ → o⋅μ and σ² → o⋅[σ²+μ²(1-o)], where o the occupancy
        @transform(:varlog  = occupancy_rescale ? :occupancy .* (:varlog .+ :meanlog.^2 .* (1 .- :occupancy)) : :varlog)
        @transform(:meanlog = occupancy_rescale ? :meanlog .* :occupancy : :meanlog)
    end
    df = innerjoin(df, sdf, on=idcolname)
    df = @chain df begin
        @subset(:varlog .> 0.)
	    @transform(:z = (:logfrequency .- :meanlog) ./ sqrt.(:varlog))
        @select(:sample_id, $(idcolname),:z)
    end
    return df
end

end # module AFD

#/ End module
