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
function compute(df::DataFrame)
    #~ Compute the (log) relative frequency of each of the "species"
    @transform!(df, :frequency = :counts ./ :nreads)
    @transform!(df, :logfrequency = log.(:frequency))
    #~ Compute some summary statistics
    nsamples = length(unique(df[!,:sample_id]))
    sdf = @chain df begin
        @by(
            :species_id,
            :occupancy = length(:sample_id) ./ nsamples,
            :meanfrequency = mean(:frequency),
            :varfrequency = var(:frequency, corrected=false),
            :meanlog = mean(:logfrequency),
            :stdlog  = std(:logfrequency, corrected=false)
        )
        @subset(:occupancy .≈ 1., :varfrequency .> 0, :stdlog .> 0)
    end
    return sdf
end

end # module Taylor
#/ End module
