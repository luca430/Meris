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
        "arxiv", "gutenberg", "rfc",
        "otu", "gowalla", "finance",
        "bcitrees", "biotime", "gtex"
    ],
    candidates = [:ParetoI, :ParetoIV, :TemperedPareto,
                  :Gamma, :LogNormal, :Weibull],
    htcandidates = [:ParetoI, :ParetoIV, :TemperedPareto]
    )
    nhtcandidates = setdiff(candidates, htcandidates)
    #~ Allocate
    summarydf = DataFrame(
        dataset=String[],
        # environment=String[],
        # fht=Float64[], K=Int[],
        nht=Int[],
        w=String,
        meanpvalue=Float64[],
        varpvalue=Float64[],
        meanexponent=Float64[],
        varexponent=Float64[],
        # pwnht=Float64[],
        # pwht=Float64[]
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
            γ = Float64[]
            nht = 0
            pv = Array{Float64}(undef, length(htrows))
            w = String[]
            #~ Compute exponents and pw
            for i in eachindex(htrows)
                pv[i] = _adf[!,:pvalue][i]
                if _adf.ntail[i] > 50 && _adf.pvalue[i] > 0.1
                    __bestmodel = _adf[!,:bestmodel][i]
                    push!(w, string(__bestmodel))
                    push!(γ, first(_fdf[!,__bestmodel][i]) + 1.0)                  
                    
                    # AICmin = _adf[!,__bestmodel][i]
                    
                    # AICs = Array(_adf[!,nhtcandidates][i,:])
                    # AICs = vcat(AICmin, AICs)
                    # lw = exp.(-0.5.*(AICs .- AICmin))
                    # push!(pwnht, lw[begin] / sum(skipmissing(lw)))

                    # AICs = Array(_adf[!,htcandidates][i,:])
                    # AICs = vcat(AICmin, AICs)
                    # lw = exp.(-0.5.*(AICs .- AICmin))
                    # push!(pwht, lw[begin] / sum(skipmissing(lw)))
                else
                    nht += 1
                end
            end

            push!(
                summarydf, [
                    "$(PROJ)-$(environment)",
                    # fht, K,
                    round(1 - nht / K, digits=2),
                    mode(w),
                    round(StatsBase.mean(pv), digits=2),
                    round(StatsBase.var(pv), digits=3),
                    round(StatsBase.mean(γ), digits=2),
                    round(StatsBase.var(γ), digits=3)
                ], promote=true
            )
        end
    end
    return summarydf
end

println(
    "lego", mean(sdf.nht), mode(sdf.w),
    mean(sdf.meanpvalue), mean(sdf.varpvalue), mean(sdf.meanexponent), mean(sdf.varexponent)
)

end # module Summarizer
#/ End module
