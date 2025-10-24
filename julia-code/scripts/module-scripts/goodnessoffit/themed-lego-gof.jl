#= Simple model to perform goodness of fit tests for the themed LEGO dataset
   For a given set of test distributions, computes the p-value for each
=#
#/ Start module
module ThemedLegoGOF

#/ Packages
using Distributions
using DataFrames, DataFramesMeta

#/ Local packages/modules
using Meris
const DATADIR = Meris.DATADIR * "macro/gof/lego/"

#################
### FUNCTIONS ###
function load_data(;
    #~ Specify variables
    minquantity = 64,         #~ Min. amount of LEGO pieces in a set
    mindistinctpieces = 32,   #~ Min. amount of distinct LEGO pieces in a set
    computefrequency = true
)
    #/ Load data
    #~ note: if `theme_id=nothing`, selects the theme with the most sets [Star Wars]
    theme_id = nothing
    legodf, themedf = Meris.LegoSampler.parse_themes(;
        minquantity=minquantity, mindistinctpieces = mindistinctpieces,
        standardize=true, returnthemes=true
    )
    #~ subselect the theme with the most sets [Star Wars `sw`]
    swlegodf = Meris.LegoSampler.select_theme(legodf, themedf)
    (computefrequency) && (@transform!(swlegodf, :frequency = :counts ./ :nreads))
    return swlegodf
end

function goodnessoffit(df::DataFrame;
    idcolname::Symbol = :component_id,
    freqcolname::Symbol = :frequency,
    minappearances::Int = 128,
    candidates = [Gamma],
    nmcsamples = 1024,
)
    #~ Compute and filter out components with insufficient no. of datapoints (appearances)
	  sdf = @chain df begin
        @by(idcolname, :nappearances = length($(idcolname)))
        @subset(:nappearances .> minappearances)
    end
    #~ Join the two dataframes to obtain only data on pieces with sufficient datapoints
    fdf = innerjoin(df, sdf, on=idcolname)
    #/ For each component, fit the distribution
    components = unique(fdf[!,:component_id])
    fitdf = DataFrame(; component_id=String[], params=[])
    for component in components
        #~ Fit using MLE
        #@TODO Make this more robust so that other methods can be used, or that distributions
        #      without a proper suff. statistic can be estimated as well.
        freq = @subset(df, :component_id .== component)[!,:frequency]
        # fitted = Distributions.fit_mle(Lomax, freq)
        fitted = Meris.MLEstimator.fit(Meris.MLEstimator.Pareto4, freq, [1.1,1.1,1.1])
        # αi, θi = params(fitted)
        push!(fitdf, [component, fitted])
        # CDF(x) = Distributions.cdf(fitted, x)
        # invCDF(x) = Distributions.quantile(fitted, x)
        # pval = Meris.GOF.estimatep(
        #     freq, CDF, invCDF, Meris.GOF.KolmogorovSmirnov, nmcsamples=nmcsamples
        # )
        # append!(pvals, pval)
    end
    return fitdf
end

end # module ThemedLegoGOF
#/ End module
