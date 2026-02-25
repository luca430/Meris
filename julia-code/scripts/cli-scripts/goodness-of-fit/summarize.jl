#= Simple script that summarizes the goodness-of-fit and candidate distributions =#
#/ Start module
module Summarizer

using StatsBase
using DataFrames, DataFramesMeta
using JLD2

using Meris:DATADIR as DATADIR
using Meris: Candies

function summarize(
    ; DIR = DATADIR*"goodness-of-fit/",
    PROJECTS = [
        "arxiv", "biotime", "finance", "otu", "rfc"
    ],
    candidates = [:GeneralizedPareto, :ParetoI, :ParetoIV, :TemperedPareto,
                  :Gamma, :LogNormal, :Weibull],
    htcandidates = [:GeneralizedPareto, :ParetoI, :ParetoIV, :TemperedPareto]
    )
    nhtcandidates = setdiff(candidates, htcandidates)
    #~ Allocate
    summarydf = DataFrame(
        dataset=String[], environment=String[],
        # fht=Float64[], K=Int[], w=Float64[],
        meanpvalue=Float64[], meanexponent=Float64[], meanpwvalue=Float64
    )
    for PROJ in PROJECTS
        #~ Load data
        jldb = JLD2.load(DIR*"$(PROJ)-candidatefits.jld2")
        fitdf = jldb["fitdf"]
        aicdf = jldb["aicdf"]
        # return fitdf, aicdf
        #~ For each `sample_id`, check which model had the lowest AIC
        @rtransform! aicdf :bestmodel = begin
            vals = collect(AsTable(candidates))
            candidates[argmin(skipmissing(vals))]
        end

        for environment in unique(fitdf.environment)
            _fdf = @subset(fitdf, :environment .== environment)
            _adf = @subset(aicdf, :environment .== environment)
            K = nrow(_fdf)

            #~ Compute fraction of models for which a heavy-tailed distribution is best
            fht = count(x -> x ∈ htcandidates, _adf.bestmodel) / nrow(_adf)
            #~ For each of the `sample_id`s where a heavy-tailed distribution is best,
            #  - extract the model parameters
            #  - compute p[w]
            htrows = findall(x -> x in htcandidates, _adf.bestmodel)
            _fdf = _fdf[htrows,:]
            _adf = _adf[htrows,:]
            
            #~ Allocate
            γ = zeros(length(htrows))
            pv = zeros(length(htrows))
            pw = zeros(length(htrows))
            #~ Compute exponents and pw
            for i in eachindex(htrows)
                γ[i] = first(_fdf[!,:ParetoI][i]) + 1.0
                AICmin = _adf[!,:ParetoI][i]
                AICs = Array(_adf[!,htcandidates][i,:])
                AICs = vcat(AICmin, AICs)
                lw = exp.(-0.5.*(AICs .- AICmin))
                pw[i] = lw[begin] / sum(skipmissing(lw))
                pv[i] = _adf[!,:pvalue][i]
            end

            push!(
                summarydf, [
                    PROJ, environment,
                    # fht, K, "ParetoIV",
                    round(StatsBase.mean(pv), digits=2),
                    round(StatsBase.mean(γ), digits=2),
                    round(StatsBase.mean(pw), digits=2)
                ], promote=true
            )
        end
    end
    return summarydf
end

end # module Summarizer
#/ End module
