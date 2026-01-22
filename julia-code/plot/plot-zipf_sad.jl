#= Module to plot Zipf, CCDF and SAD considering a filter to remove low counts =#
#/ Start module
module SADPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta
using StatsBase, JLD2

#/ Modules
import Meris

#################
### FUNCTIONS ###
function plot_taylor(;
        ZIPFDIR=Meris.DATADIR * "macro/zipf/",
        relative_counts=false,
        savefig=false,
        figname="zipf.png"
    )
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    count_label = relative_counts ? "ν" : "n"

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(2.25*width,height), figure_padding=(2,4,2,14))
    
    ax1 = Axis(
        fig[1,1],
        xlabel=L"\text{rank}", xlabelsize=11,
        ylabel=L"\log_{10}(\text{%$count_label})", ylabelsize=11,
        yscale=log10, xscale=log10
    )
    
    ax2 = Axis(
        fig[1,2],
        xlabel=L"\log_{10}(\text{%$count_label})", xlabelsize=11,
        ylabel=L"\text{counts pdf}", ylabelsize=11,
        yscale=log10
    )

    ax3 = Axis(
        fig[1,3],
        xlabel=L"\alpha", xlabelsize=11,
        ylabel=L"\alpha \text{ pdf}", ylabelsize=11,
        yscale=log10
    )

    # Load data
    rfc = JLD2.load(ZIPFDIR * "rfc.jld2")["figure"]
    gtex = JLD2.load(ZIPFDIR * "gtex.jld2")["figure"]
    otu = JLD2.load(ZIPFDIR * "otu.jld2")["figure"]
    otu2 = JLD2.load(ZIPFDIR * "otu2.jld2")["figure"]
    otu3 = JLD2.load(ZIPFDIR * "otu3.jld2")["figure"]

    # Plot ax1
    # scatter!(ax1, rfc.ax1.scatterx, rfc.ax1.scattery, color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4)
    # lines!(ax1, rfc.ax1.linex, rfc.ax1.liney,
    #         linewidth=0.8, linestyle=:dash, color=:black, label=L"\xi = %$(round(1 /rfc.α, digits=2))")

    # scatter!(ax1, gtex.ax1.scatterx, gtex.ax1.scattery, color=:white, strokecolor=colors[2], markersize=4, strokewidth=0.4)
    # lines!(ax1, gtex.ax1.linex, gtex.ax1.liney,
    #         linewidth=0.8, linestyle=:dash, color=:black, label=L"\xi = %$(round(1 /gtex.α, digits=2))")

    scatter!(ax1, otu.ax1.scatterx, otu.ax1.scattery, color=:white, strokecolor=colors[3], markersize=4, strokewidth=0.4)
    lines!(ax1, otu.ax1.linex, otu.ax1.liney,
            linewidth=0.8, linestyle=:dash, color=colors[3], label=L"\xi = %$(round(1 /otu.α, digits=2))")

    scatter!(ax1, otu2.ax1.scatterx, otu2.ax1.scattery, color=:white, strokecolor=colors[4], markersize=4, strokewidth=0.4)
    lines!(ax1, otu2.ax1.linex, otu2.ax1.liney,
            linewidth=0.8, linestyle=:dash, color=colors[4], label=L"\xi = %$(round(1 /otu2.α, digits=2))")

    scatter!(ax1, otu3.ax1.scatterx, otu3.ax1.scattery, color=:white, strokecolor=colors[5], markersize=4, strokewidth=0.4)
    lines!(ax1, otu3.ax1.linex, otu3.ax1.liney,
            linewidth=0.8, linestyle=:dash, color=colors[5], label=L"\xi = %$(round(1 /otu3.α, digits=2))")
    axislegend(ax1, position=:lb)

    # Plot ax2
    # scatter!(ax2, rfc.ax2.scatterx, rfc.ax2.scattery, color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4)
    # lines!(ax2, rfc.ax2.linex, rfc.ax2.liney,
    #         linewidth=0.8, linestyle=:dash,  color=:black, label=L"\gamma = %$(round(rfc.α + 1, digits=2))")

    # scatter!(ax2, gtex.ax2.scatterx, gtex.ax2.scattery, color=:white, strokecolor=colors[2], markersize=4, strokewidth=0.4)
    # lines!(ax2, gtex.ax2.linex, gtex.ax2.liney,
    #         linewidth=0.8, linestyle=:dash,  color=:black, label=L"\gamma = %$(round(gtex.α + 1, digits=2))")
    
    scatter!(ax2, otu.ax2.scatterx, otu.ax2.scattery, color=:white, strokecolor=colors[3], markersize=4, strokewidth=0.4)
    lines!(ax2, otu.ax2.linex, otu.ax2.liney,
            linewidth=0.8, linestyle=:dash,  color=colors[3], label=L"\gamma = %$(round(otu.α + 1, digits=2))")

    scatter!(ax2, otu2.ax2.scatterx, otu2.ax2.scattery, color=:white, strokecolor=colors[4], markersize=4, strokewidth=0.4)
    lines!(ax2, otu2.ax2.linex, otu2.ax2.liney,
            linewidth=0.8, linestyle=:dash,  color=colors[4], label=L"\gamma = %$(round(otu2.α + 1, digits=2))")

    scatter!(ax2, otu3.ax2.scatterx, otu3.ax2.scattery, color=:white, strokecolor=colors[5], markersize=4, strokewidth=0.4)
    lines!(ax2, otu3.ax2.linex, otu3.ax2.liney,
            linewidth=0.8, linestyle=:dash,  color=colors[5], label=L"\gamma = %$(round(otu3.α + 1, digits=2))")
    axislegend(ax2, position=:lb)

    # Plot ax3
    # scatter!(ax3, rfc.ax3.scatterx, rfc.ax3.scattery ./ maximum(rfc.ax3.scattery), color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4)
    # scatter!(ax3, gtex.ax3.scatterx, gtex.ax3.scattery ./ maximum(gtex.ax3.scattery), color=:white, strokecolor=colors[2], markersize=4, strokewidth=0.4)
    scatter!(ax3, otu.ax3.scatterx, otu.ax3.scattery ./ maximum(otu.ax3.scattery), color=:white, strokecolor=colors[3], markersize=4, strokewidth=0.4)
    scatter!(ax3, otu2.ax3.scatterx, otu2.ax3.scattery ./ maximum(otu2.ax3.scattery), color=:white, strokecolor=colors[4], markersize=4, strokewidth=0.4)
    scatter!(ax3, otu3.ax3.scatterx, otu3.ax3.scattery ./ maximum(otu3.ax3.scattery), color=:white, strokecolor=colors[5], markersize=4, strokewidth=0.4)

    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end






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
        xlabel=L"\text{rank}", xlabelsize=11,
        ylabel=L"\log_{10}(\text{%$count_label})", ylabelsize=11,
        xscale=log10
    )
    
    ax2 = Axis(
        fig[1,2],
        xlabel=L"\log_{10}(\text{%$count_label})", xlabelsize=11,
        ylabel=L"\text{pdf}", ylabelsize=11,
        yscale=log10
    )

    if !aggregate
        cdf = deepcopy(df)
        (relative_counts) && (cdf.counts ./= cdf.nreads)
        zipf_df = zipf(cdf)
        
        samples = unique(df.sample_id)
        ε_vec = []
        α_vec = []
        max_count_vec = []
        max_rank_vec = []
        for (i,sample) in enumerate(samples[samples_idx])
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
            
            if filter
                mask = counts .> ε
                ranks = ranks[mask]
                counts = counts[mask]
            end

            push!(max_count_vec, maximum(counts))
            push!(max_rank_vec, maximum(ranks))
            
            scatter!(
                ax1, ranks, log10.(counts),
                color=:white, strokecolor=colors[i], markersize=4, strokewidth=0.4
            )
        
            x, pdf = Meris.DataTools.make_hist(log10.(counts); nbins=nbins)
            scatter!(
                ax2, x, pdf ./ (α * ε ^ α),
                color=:white, strokecolor=colors[i], markersize=4, strokewidth=0.4
            )
        end

        pareto = Meris.ParetoDistribution.ParetoI(mean(α_vec), mean(ε_vec))
        paretox = relative_counts ? collect(mean(ε_vec):1e-4:minimum(max_count_vec)) : collect(mean(ε_vec):minimum(max_count_vec))
        Z = sum(count(x -> x > minimum(ε_vec), paretox)) / length(paretox)
        lines!(ax1, mean(max_rank_vec) .* Meris.ParetoDistribution.ccdf.(pareto, paretox) .* Z, log10.(paretox),
            linewidth=1, linestyle=:dash, color=:black, label=L"\alpha = %$(round(1 /mean(α_vec), digits=2))")
        axislegend(ax1, position=:lb)
        lines!(ax2, log10.(paretox), paretox .^ (-mean(α_vec)) .* log(10),
            linewidth=0.8, linestyle=:dash,  color=:black, label=L"\gamma = %$(round(mean(α_vec) + 1, digits=2))")
        axislegend(ax2, position=:lb)
        
    else
        cdf = deepcopy(df)
        samples = unique(df.sample_id)
        filter!(row -> row.sample_id in samples[samples_idx], cdf)
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
        
        scatter!(
            ax1, ranks, log10.(counts),
            color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4
        )
    
        x, pdf = Meris.DataTools.make_hist(log10.(counts); nbins=nbins)
        scatter!(
            ax2, x, pdf ./ (α * ε ^ α),
            color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4
        )
        
        pareto = Meris.ParetoDistribution.ParetoI(α, ε)
        paretox = relative_counts ? collect(ε:1e-4:maximum(counts)) : collect(ε:maximum(counts))
        lines!(ax1, maximum(ranks) .* Meris.ParetoDistribution.ccdf.(pareto, paretox), log10.(paretox),
            linewidth=1, linestyle=:dash,  color=:black, label=L"\gamma = %$(round((1 / α), digits=2))")
        axislegend(ax1, position=:lb)
        lines!(ax2, log10.(paretox), paretox .^ (-α) .* log(10),
            linewidth=0.8, linestyle=:dash,  color=:black, label=L"\alpha = %$(round(α + 1, digits=2))")
        axislegend(ax2, position=:lb)
        
    end

    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))

    return fig
end



end # module SADPlotter
#/ End module
