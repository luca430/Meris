#= Module to plot Zipf, CCDF and SAD considering a filter to remove low counts =#
#/ Start module
module SADPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta
using StatsBase, Random

#/ Modules
import Meris

#################
### FUNCTIONS ###
function plot(
        df;
        samples_idx=nothing,
        relative_counts=false,
        aggregate=false,
        filter=true,
        xmins=nothing,
        savefig = false,
        figname = true,
        nbins=30
    )

    Random.seed!(1234)
    default_idx = rand(1:length(unique(df.sample_id)), 10)
    samples_idx = isnothing(samples_idx) ? default_idx : samples_idx
    count_label = relative_counts ? "ν" : "n"
    
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(1.5*width,height), figure_padding=(2,4,2,14))
    
    ax1 = Axis(
        fig[1,1],
        xlabel=L"\text{rank}",
        ylabel=L"\log_{10}(\text{%$count_label})",
        xscale=log10
    )
    
    ax2 = Axis(
        fig[1,2],
        xlabel=L"\log_{10}(\text{%$count_label})",
        ylabel=L"\text{pdf}",
        yscale=log10
    )

    if !aggregate
        cdf = deepcopy(df)
        (relative_counts) && (cdf.counts ./= cdf.nreads)
        zipf_df = zipf(cdf)
        
        samples = unique(df.sample_id)
        ε_vec = []
        α_vec = []
        Z_vec = []
        max_count_vec = []
        max_rank_vec = []
        for sample in samples[samples_idx]
            sdf = zipf_df[zipf_df.sample_id .== sample, :]
            ranks, counts = sdf.ranks, sdf.counts
            p = sortperm(ranks)          # permutation that sorts ranks ascending
            ranks  = ranks[p]
            counts = counts[p]

            fit = Meris.Powerlaw.fitPareto(sdf.counts; xmins=xmins, minsamples=50)
            ε = fit.Pareto.ε
            push!(ε_vec, ε)
            α = fit.Pareto.α
            push!(α_vec, α)
            Z = sum(count(x -> x > ε, counts)) / length(counts)
            push!(Z_vec, Z)

            if filter
                mask = counts .> ε
                ranks = ranks[mask]
                counts = counts[mask]
            end

            push!(max_count_vec, maximum(counts))
            push!(max_rank_vec, maximum(ranks))
            
            scatter!(ax1, ranks, log10.(counts), markersize=1)
        
            x, pdf = Meris.DataTools.make_hist(log10.(counts); nbins=nbins)
            scatter!(ax2, x, pdf, markersize=1)
        end

        pareto = Meris.ParetoDistribution.ParetoI(mean(α_vec), minimum(ε_vec))
        paretox = relative_counts ? collect(minimum(ε_vec):1e-4:minimum(max_count_vec)) : collect(minimum(ε_vec):minimum(max_count_vec))
        lines!(ax1, mean(max_rank_vec) .* Meris.ParetoDistribution.ccdf.(pareto, paretox), log10.(paretox),
            linewidth=1, linestyle=:dash, color=:black, label=L"\gamma = %$(round(1 /mean(α_vec), digits=2))")
        axislegend(ax1)
        lines!(ax2, log10.(paretox), log(10) .* paretox .* Meris.ParetoDistribution.pdf.(pareto, paretox),
            linewidth=0.8, linestyle=:dash,  color=:black, label=L"\alpha = %$(round(mean(α_vec) + 1, digits=2))")
        axislegend(ax2)
        
    else
        cdf = deepcopy(df)
        (relative_counts) && (cdf.counts ./= cdf.nreads)
        zipf_df = agg_zipf(cdf)
        
        ranks, counts = zipf_df.ranks, zipf_df.counts
        p = sortperm(ranks)          # permutation that sorts ranks ascending
        ranks  = ranks[p]
        counts = counts[p]

        fit = Meris.Powerlaw.fitPareto(zipf_df.counts; xmins=xmins, minsamples=50)
        ε = fit.Pareto.ε
        α = fit.Pareto.α
        Z = sum(count(x -> x > ε, counts)) / length(counts)

        if filter
            mask = counts .> ε
            ranks = ranks[mask]
            counts = counts[mask]
        end
        
        scatter!(ax1, ranks, log10.(counts), markersize=1)
    
        x, pdf = Meris.DataTools.make_hist(log10.(counts); nbins=nbins)
        scatter!(ax2, x, pdf, markersize=1)
        
        pareto = Meris.ParetoDistribution.ParetoI(α, ε)
        paretox = relative_counts ? collect(ε:1e-4:maximum(counts)) : collect(ε:maximum(counts))
        lines!(ax1, maximum(ranks) .* Meris.ParetoDistribution.ccdf.(pareto, paretox), log10.(paretox),
            linewidth=1, linestyle=:dash,  color=:black, label=L"\gamma = %$(round((1 / α), digits=2))")
        axislegend(ax1)
        lines!(ax2, log10.(paretox), log(10) .* paretox .* Meris.ParetoDistribution.pdf.(pareto, paretox),
            linewidth=0.8, linestyle=:dash,  color=:black, label=L"\alpha = %$(round(α + 1, digits=2))")
        axislegend(ax2)
        
    end

    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))

    return fig
end


########################
### HELPER FUNCTIONS ###
function zipf(df)
    return @chain df begin
    @groupby(:sample_id)
    @combine(:ranks = tiedrank(-:counts), :counts, :nreads)
    end
end

function agg_zipf(df)
    return @chain df begin
        @groupby(:component_id)
        @combine(:sample_id, :agg_counts = sum(:counts), :agg_nreads = sum(:nreads))
        @transform(
            :ranks = tiedrank(-:agg_counts),
            :counts = :agg_counts,
            :nreads = :agg_nreads
        )
    end
end

##############################

end # module SADPlotter
#/ End module
