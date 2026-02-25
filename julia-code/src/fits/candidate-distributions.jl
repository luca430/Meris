#= Small module that contains candidate distributions and appropriate methods =#
#/ Start module
module Candies

using Distributions
using Meris: MDistributions

#################
### FUNCTIONS ###
"""
Simple function that returns dictionary of candidates and appropriate method calls that
can be iterated over
"""
function getcandidates(
    ; candidates=[:GeneralizedPareto, :ParetoI, :ParetoIV, :TemperedPareto,
                  :Gamma, :LogNormal, :Weibull]
    )
    #~ List all available candidates here
    allcandidates = Dict(
        :GeneralizedPareto => (
            ; f = MDistributions.GeneralizedPareto,
            fit = (f, data, εs) -> MDistributions.fit(f, data, εs),
            p = (f) -> MDistributions.params(f),
            logpdf = (f, x) -> MDistributions.logpdf.(f, x),
            computepvalue = MDistributions.computepvalue,
            dataframeentry = Tuple{Float64,Float64,Float64}
        ),
        :ParetoI => (
            ; f = MDistributions.ParetoI,
            fit = (f, data, εs) -> MDistributions.fit(f, data, εs),
            p = (f) -> MDistributions.params(f),
            logpdf = (f, x) -> log.(MDistributions.pdf.(f, x)),
            computepvalue = MDistributions.computepvalue,
            dataframeentry = Tuple{Float64,Float64}
        ),
        :ParetoIV => (
            ; f = MDistributions.ParetoIV,
            fit = (f, data, εs) -> MDistributions.fit(f, data, 1.0, εs),
            p = (f) -> MDistributions.params(f),
            logpdf = (f, x) -> MDistributions.logpdf.(f, x),
            #~ Defaults to Burr/Lomax distribution with β=1
            computepvalue = (P, data, εs; β=1.0, weighted=false, rng=rng) ->
                MDistributions.computepvalue(P, data, β, εs; weighted=weighted),
            dataframeentry = Tuple{Float64,Float64,Float64,Float64}
        ),
        :TemperedPareto => (
            ; f = MDistributions.TemperedPareto,
            fit = (f, data, εs) -> MDistributions.fit(f, data, εs),
            p = (f) -> MDistributions.params(f),
            logpdf = (f, x) -> MDistributions.logpdf.(f, x),
            computepvalue = MDistributions.computepvalue,
            dataframeentry = Tuple{Float64,Float64,Float64}
        ),
        :Gamma => (
            ; f = Distributions.Gamma,
            fit = (f, data, εs) -> Distributions.fit_mle(f, data),
            p = (f) -> Distributions.params(f),
            logpdf = (f, x) -> log.(Distributions.pdf.(f, x)),
            dataframeentry = Tuple{Float64,Float64}
        ),
        :LogNormal => (
            ; f = Distributions.LogNormal,
            fit = (f, data, εs) -> Distributions.fit_mle(f, data),
            p = (f) -> Distributions.params(f),
            logpdf = (f, x) -> log.(Distributions.pdf.(f, x)),
            dataframeentry = Tuple{Float64,Float64}
        ),
        :Weibull => (
            ; f = Distributions.Weibull,
            fit = (f, data, εs) -> Distributions.fit_mle(f, data),
            p = (f) -> Distributions.params(f),
            logpdf = (f, x) -> log.(Distributions.pdf.(f, x)),
            dataframeentry = Tuple{Float64,Float64}
        )
    )
    #~ Filter based on desired candidates
    return filter(entry -> entry.first in candidates, allcandidates)
end

end # module CDistributions
#/ End module
