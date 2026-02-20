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

