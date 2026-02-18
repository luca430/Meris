#= Simple module with function to fit candidate distributions to data
   Here, candidate distributions may be listed, including their methods to call and fit them.
   Subsequently, a dataframe can be initialized with the appropriate fields.
   This is done manually as many candidate distributions have a different set of parameters
   and specified function to call and fit.
=#

#/ Start module
module Candidates

using DataFrames
using Distributions

using Meris: MDistributions

#= Short-list of candidate distributions =#
function define_candidates()
    candidates = Dict(
        :ParetoI => (
            ; f = MDistributions.ParetoI,
            fit = (f, data, εs) -> MDistributions.fit(f, data; εs=εs),
            p = (f) -> MDistributions.params(f)
        ),
        :ParetoIV => (
            ; f = MDistributions.ParetoIV,
            fit = (f, data, εs) -> MDistributions.fit(f, data; εs=εs),
            p = (f) -> MDistributions.params(f)
        ),
        :TemperedPareto => (
            ; f = MDistributions.TemperedPareto,
            fit = (f, data, εs) -> MDistributions.fit(f, data; εs=εs),
            p = (f) -> MDistributions.params(f)
        ),
        :Gamma => (
            ; f = Distributions.Gamma,
            fit = (f, data, εs) -> Distributions.fit_mle(f, data),
            p = (f) -> Distributions.params(f)
        ),
        :LogNormal => (
            ; f = Distributions.LogNormal,
            fit = (f, data, εs) -> Distributions.fit_mle(f, data),
            p = (f) -> Distributions.params(f)
        ),
        :Weibull => (
            ; f = Distributions.Weibull,
            fit = (f, data, εs) -> Distributions.fit_mle(f, data),
            p = (f) -> Distributions.params(f)
        )
    )
    return candidates
end

function initialize_dataframe()
    fitdf = DataFrame(
        environment=String[], sample_id=String[],
        ParetoI=Tuple{Float64,Float64},
        ParetoIV=Tuple{Float64,Float64,Float64,Float64},
        TemperedPareto=Tuple{Float64,Float64,Float64},
        Gamma=Tuple{Float64,Float64},
        LogNormal=Tuple{Float64,Float64},
        Weibull=Tuple{Float64,Float64}
    )
    return fitdf
end

function fit_candidates(data::DataFrame; nε::Int = 64)
    #~ Get candidates and allocate
    candidates = define_candidates()
    fitdf = initialize_dataframe()
    
    nsamples = length(unique(data.sample_id))
    n = 0
    for sampledf in groupby(data, :sample_id)
        n += 1
        __id = first(sampledf.sample_id)
        println("Fitting candidate distributions for $(__id) [$(n)/$(nsamples)]...")
        frequencies = collect(sampledf.frequency)
        #~ Compute admissible ε for distributions from the Pareto family
        νs = unique(sort(frequencies))
        εs = range(νs[2], νs[end], nε) |> collect
        #~ Fit candidate distributions
        fits = Dict{Symbol,Any}()
        for (name, distribution) in candidates            
            fits[name] = distribution.p(distribution.fit(distribution.f, frequencies, εs))
        end
        push!(
            fitdf, merge((environment = first(sampledf.domain), sample_id = __id), fits),
            promote=true
        )
    end
    return fitdf
end

end # module Candidates
#/ End module

