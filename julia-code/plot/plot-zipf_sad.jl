#= Module to plot Zipf, CCDF and SAD considering a filter to remove low counts =#
#/ Start module
module SADPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta
using StatsBase, JLD2
using Colors, ColorTypes

#/ Modules
import Meris
import Meris.MDistributions as MDist

#################
### FUNCTIONS ###
function plot(;
        color_num=1,
        ZIPFDIR=Meris.DATADIR * "macro/zipf/",
        relative_counts=false,
        savefig=false,
        figname="zipf.png"
    )
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    base = MakiePublication.COLORS[1][color_num]
    base_hsl = HSL(base)  # convert to HSL
    colors = [HSL(base_hsl.h, base_hsl.s, l) for l in range(0.1, 0.8, length=5)]

    count_label = relative_counts ? "ν" : "n"

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(2.25*width,height), figure_padding=(2,4,2,14))

    ax1 = Axis(
        fig[1,1],
        xlabelsize=11, xticklabelsize=6,
        ylabel=L"\text{CAD exponent } \gamma", ylabelsize=11
    )

    ax2 = Axis(
        fig[1,2],
        xlabel=L"\text{Sample size } N", xlabelsize=11,
        ylabel=L"\text{Vocabulary size } V", ylabelsize=11,
        xscale=log10, yscale=log10
    )
    
    ax3 = Axis(
        fig[1,3],
        xlabel=L"\text{rank } r", xlabelsize=11,
        ylabel=L"\log_{10}(\text{%$count_label})", ylabelsize=11,
        yscale=log10, xscale=log10
    )

    # Load data
    otu = JLD2.load(ZIPFDIR * "otu.jld2")["out"]

    # Axis 1: Boxplot
    labels = collect(keys(otu.boxplot))
    vals = collect(values(otu.boxplot))
    
    for (i, v) in enumerate(vals)
        boxplot!(ax1, fill(i, length(v)), v, color=(colors[i], 0.8), markersize=4, whiskerlinewidth=0.8, medianlinewidth=0.8)
    end
    ax1.xticks = (1:length(labels), labels)

    # Axis 2: Heaps' law
    for (i,v) in enumerate(values(otu.heaps))
        Δ = Int(ceil(length(v.samplesize) / 50))
        scatter!(ax2, v.samplesize[1:Δ:end], v.vocabsize[1:Δ:end], color=:white, strokecolor=colors[i], markersize=5, strokewidth=0.5)
    end

    # Axis 3: Zipf's law
    for (i,(v,f)) in enumerate(zip(values(otu.zipf), values(otu.fit)))
        Δ = Int(ceil(length(v.ranks) / 60))
        scatter!(ax3, v.ranks[1:Δ:end], v.counts[1:Δ:end], color=:white, strokecolor=colors[i], markersize=5, strokewidth=0.5)
        
        func = MDist.ParetoI(f.α, f.ε)
        x = 10 .^ collect(log10(minimum(v.counts)):1e-2:log10(maximum(v.counts)))
        lines!(ax3, MDist.ccdf.(func, x) .* maximum(v.ranks), x, linestyle=:dash, color=:black, linewidth=0.8)
    end

    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

### HELPER ###

function shades(color::Colorant, n; lmin=0.35, lmax=0.8)
    hsl = HSL(color)
    ls = range(lmin, lmax; length=n)
    return [RGB(HSL(hsl.h, hsl.s, l)) for l in ls]
end

end # module SADPlotter
#/ End module
