#= Simple module to compute the mean and frequency for a given DataFrame
   Useful for investigating Taylor's law
=#
#/ Start module
module Taylor

#/ Packages
using DataFrames, DataFramesMeta
using StatsBase
using LsqFit

using Meris

#################
### FUNCTIONS ###
""
function compute(df::DataFrame, idcolname; minoccupancy::Float64=0., maxfrequency::Float64=1.)
    #~ Compute the (log) relative frequency of each of the "species"
    @transform!(df, :frequency = :counts ./ :nreads)
    @transform!(df, :logfrequency = log.(:frequency))
    #~ Compute some summary statistics
    nsamples = length(unique(df[!, :sample_id]))
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
        @subset(:occupancy .>= minoccupancy, :s .> 0.)
        #~ Select only those with a rel. frequency below the given maximum, and positive variance
        #  note: The first filters out those where, roughly, a Poisson distribution approximates
        #        a Binomial distribution well.
        @subset(:m .< maxfrequency)
        #~ Take the occupation number into account
        #~ this means that μ → o⋅μ and σ² → o⋅[σ²+μ²(1-o)], where o the occupancy
        #! note that here σ²=s
        @transform(:omeanfrequency = :m .* :occupancy)
        @transform(:ovarfrequency = :occupancy .* (:s .+ :m .^ 2 .* (1 .- :occupancy)))
        #~ Perform a log-transform on the mean-frequency (needed for lognormal)
        @rename(:meanfrequency = :m, :varfrequency = :s)
        @select(
            ($idcolname),
            :meanfrequency, :varfrequency,
            :omeanfrequency, :ovarfrequency,
            :noccurrences, :occupancy,
            :varm, :vars,
            :errorcov, :errorcorr
        )
    end
    return sdf
end

### Alternative function to compute TL from a count matrix
### Note: it requires to first convert a df into a count matrix
function compute2(counts; binned=true, nbins=30)

    means = vcat(mean(counts, dims=1)...)
    vars = vcat(var(counts, dims=1)...)

    mask = vars .> 0
    log_means = log10.(means[mask])
    log_vars = log10.(vars[mask])

    if binned
        log_means, log_vars = binned_average(log_means, log_vars; nbins=nbins)
        log_means, log_vars = clean_log(log_means, log_vars)
    end

    # model(x, p) = p[1] .+ p[2] .* x
    # fit = curve_fit(model, log_means, log_vars, [0.0, 2.0])
    # weights = ones(length(log_means))
    # fit = Meris.StraightLine.weightedyorkfit(log_means, log_vars, weights, weights)

    return (log_means, log_vars)
end


# Helper functions
function binned_average(x, y; nbins=20)
    edges = range(minimum(x), stop=maximum(x), length=nbins + 1)
    bin_indices = searchsortedfirst.(Ref(edges), x) .- 1  # get bin index for each x
    bin_indices = clamp.(bin_indices, 1, nbins)  # ensure indices are within range

    keep = [any(bin_indices .== i) for i in 1:nbins]
    y_mean = [mean(y[bin_indices .== i]) for i in 1:nbins if keep[i]]
    x_center = [(edges[i] + edges[i+1]) / 2 for i in 1:nbins if keep[i]]

    return x_center, y_mean
end

function clean_log(x, y)
    mask = isfinite.(x) .& isfinite.(y)
    return x[mask], y[mask]
end

end # module Taylor
#/ End module
