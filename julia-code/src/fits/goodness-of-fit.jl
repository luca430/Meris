#= Module for goodness of fit tests =#
#/ Start module
module OhMyGoodness

#/ Packages
using Distributions
using DataFrames, DataFramesMeta

#~ Local Meris modules
using Meris: Candies   #~ candidate distributions

#################
### FUNCTIONS ###
"""
Fit a set of candidate distributions to the data
"""
function fit_candidates(
    data::DataFrame, classcolname::Symbol;
    candidates = [:GeneralizedPareto, :ParetoI, :ParetoIV, :TemperedPareto,
                  :Gamma, :LogNormal, :Weibull],
    nε::Int = 100,
    preject = 0.1
    )
    candidates = Candies.getcandidates(; candidates=candidates)
    fitdf = initialize_fitdataframe(candidates)
    aicdf = initialize_aicdataframe(candidates)

    #~ For each `sample_id` in the data DataFrame, fit all candidate distributions
    nsamples = length(unique(data.sample_id))
    n = 0
    nrejects = 0
    for sampledf in groupby(data, :sample_id)
        p = missing
        xmin = nothing
        ntail = nothing
        n += 1
        __id = first(sampledf.sample_id)
        frequencies = collect(sampledf.frequency)
        #~ Compute [log-spaced] admissible ε for distributions from the Pareto family
        νs = log.(unique(sort(frequencies)))
        (length(νs) < 3) && (continue)
        εs = exp10.(range(νs[2], νs[end], nε) |> collect)
        #~ Establish heavy-tail by fitting generalized Pareto distribution
        try
            ht = candidates[:ParetoI]
            paretofit = ht.fit(ht.f, frequencies, εs)
            #~ Compute p value
            p = ht.computepvalue(paretofit, frequencies, εs)
            if p < preject
                nrejects += 1
                continue
            end            
            #~ Filter data w.r.t. ε onwards, otherwise (i)
            #  - the log-likelihood blows up
            #  - the comparison is unfair as candidates have different domains
            frequencies = frequencies[frequencies .>= paretofit.ε]
            xmin = paretofit.ε
            ntail = length(frequencies)
        catch DomainError
            # Sometimes fitting a Pareto-like distribution is troublesome, so catch
            # any potential errors here and simply skip the source.
            continue
        end
        #~ Once we are here, a heavy-tailed distribution is _not_ rejected, so now we
        #  compare it with the other remaining candidate distributions.
        #  So from this point on, we filter frequencies by ε of the test distribution
        #~ Fit candidate distributions
        print("Fitting candidate distributions for sample $(__id) [$(n)/$(nsamples)]...\r")
        fits = Dict{Symbol,Any}()
        AICs = Dict{Symbol,Any}()
        for (name, distribution) in candidates
            try
                __fit = distribution.fit(distribution.f, frequencies, xmin)
                fits[name] = distribution.p(__fit)
                #~ Compute Akaike information criterion
                L = -sum(distribution.logpdf(__fit, frequencies))
                u = length(distribution.p(__fit))
                AICs[name] = 2*u + 2*L
            catch
                fits[name] = missing
                AICs[name] = missing
            end
        end
        fitmrg = merge(
            (environment = first(sampledf[!,classcolname]), sample_id = __id), fits
        )
        push!(fitdf, fitmrg, promote=true)
        aicmrg = merge(
            (environment = first(sampledf[!,classcolname]),
             sample_id = __id, pvalue=p, ntail=ntail), AICs
        )
        push!(aicdf, aicmrg, promote=true)
    end
    println("\nDone.")
    return fitdf, aicdf
end

########################
### HELPER FUNCTIONS ###
"Initialize DataFrame for fitting purposes"
function initialize_fitdataframe(candidates::Dict)
    df = DataFrame(environment=String[], sample_id=String[])
    for (candidate, properties) in candidates
        df[:, candidate] = Vector{properties.dataframeentry}(undef, nrow(df))
    end
    return df
end

function initialize_aicdataframe(candidates::Dict)
    df = DataFrame(environment=String[], sample_id=String[], pvalue=Float64[], ntail=Int[])
    for (candidate, _) in candidates
        df[:, candidate] = Vector{Float64}(undef, nrow(df))
    end
    return df
end

end # module Goodness
#/ End module
