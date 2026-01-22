module Zipf

using Meris
using DataFrames, DataFramesMeta, StatsBase
using JLD2, Random

function compute(
        df;
        relative_counts=false,
        filter=true,
        xmins=nothing,
        nbins=30,
        save=true,
        filename="zipf.jld2"
    )

    Random.seed!(1234)

    agg_df = aggregate_samples(df)
    (relative_counts) && (agg_df.counts ./= agg_df.nreads)

    α_vec, ε_vec = fit_samples(agg_df; xmins=xmins)
    α_x, α_pdf = Meris.DataTools.make_hist(α_vec; nbins=40)

    agg_df.samples_id .= "samp"
    α, ε = fit_samples(agg_df; xmins=xmins)
    α, ε = α[1], ε[1]

    agg_df.ranks .= tiedrank(-agg_df.counts)
    ranks, counts = agg_df.ranks, agg_df.counts
    p = sortperm(ranks)          # permutation that sorts ranks ascending
    ranks  = ranks[p]
    counts = counts[p]
    
    if filter
        mask = counts .> ε
        ranks = ranks[mask]
        counts = counts[mask]
    end

    x, pdf = Meris.DataTools.make_hist(log10.(counts); nbins=nbins)
    pdf = pdf ./ (α * ε ^ α)

    pareto = Meris.ParetoDistribution.ParetoI(α, ε)
    paretox = relative_counts ? collect(ε:1e-4:maximum(counts)) : collect(ε:maximum(counts))

    ax1 = (
        scatterx = ranks,
        scattery = counts,
        linex = maximum(ranks) .* Meris.ParetoDistribution.ccdf.(pareto, paretox),
        liney = paretox
        )
    ax2 = (
        scatterx = x,
        scattery = pdf,
        linex = log10.(paretox),
        liney = paretox .^ (-α) .* log(10)
        )
    ax3 = (
        scatterx = α_x,
        scattery = α_pdf
    )

    figure = (ax1 = ax1, ax2 = ax2, ax3 = ax3, α = α)
    (save) && (@save filename figure)
    
    return figure
end

#### HELPER ####

function aggregate_samples(df)
    return @chain df begin
        @groupby(:component_id)
        @combine(:sample_id, :agg_counts = sum(:counts), :agg_nreads = sum(:nreads))
        @transform(
            :counts = :agg_counts,
            :nreads = :agg_nreads
        )
    end
end
    
function fit_samples(df; samples_idx=nothing, xmins=nothing)
    
    samples = unique(df.sample_id)
    samples_idx = isnothing(samples_idx) ? collect(1:length(samples)) : samples_idx

    ε_vec = []
    α_vec = []
    for sample in samples[samples_idx]
        sdf = df[df.sample_id .== sample, :]

        fit = Meris.Powerlaw.fitPareto(sdf.counts; xmins=xmins, minsamples=50)
        ε = fit.Pareto.ε
        push!(ε_vec, ε)
        α = fit.Pareto.α
        push!(α_vec, α)
    end
    
    return (α_vec, ε_vec)
end

end # End module