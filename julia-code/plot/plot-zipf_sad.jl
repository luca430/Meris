#= Module to plot Zipf, CCDF and SAD considering a filter to remove low counts =#
#/ Start module
module SADPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta
using StatsBase, JLD2
using Colors

#/ Modules
import Meris

#################
### FUNCTIONS ###
function plot_taylor(;
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
    colors = shades(base, 5)

    count_label = relative_counts ? "ν" : "n"

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(2.25*width,height), figure_padding=(2,4,2,14))

    ax1 = Axis(
        fig[1,1],
        xlabel=L"\log_{10}(\text{%$count_label})", xlabelsize=11,
        ylabel=L"\text{counts pdf}", ylabelsize=11,
        yscale=log10,
        limits=(-5,-0.9,5e0,1e5)
    )
    
    ax2 = Axis(
        fig[1,2],
        xlabel=L"\text{rank}", xlabelsize=11,
        ylabel=L"\log_{10}(\text{%$count_label})", ylabelsize=11,
        yscale=log10, xscale=log10,
        limits=(5e-4,1,1,2e2)
    )

    ax3 = Axis(
        fig[1,3],
        xlabel=L"\text{Sample size } N", xlabelsize=11,
        ylabel=L"\text{Vocabulary size } V", ylabelsize=11
    )

    # Load data
    otu = JLD2.load(ZIPFDIR * "otu.jld2")["figure"]
    otu2 = JLD2.load(ZIPFDIR * "otu2.jld2")["figure"]
    otu3 = JLD2.load(ZIPFDIR * "otu3.jld2")["figure"]
    otu5 = JLD2.load(ZIPFDIR * "otu5.jld2")["figure"]

    # Plot ax1
    scatter!(ax1, otu.ax1.scatterx, otu.ax1.scattery, color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4)
    scatter!(ax1, otu2.ax1.scatterx, otu2.ax1.scattery, color=:white, strokecolor=colors[2], markersize=4, strokewidth=0.4)
    scatter!(ax1, otu3.ax1.scatterx, otu3.ax1.scattery, color=:white, strokecolor=colors[3], markersize=4, strokewidth=0.4)
    scatter!(ax1, otu5.ax1.scatterx, otu5.ax1.scattery, color=:white, strokecolor=colors[4], markersize=4, strokewidth=0.4)
    lines!(ax1, otu5.ax1.linex, otu5.ax1.liney,
            linewidth=1, linestyle=:dash,  color=:black)

    # Plot ax2
    scatter!(ax2, otu.ax2.scatterx, otu.ax2.scattery, color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4)
    scatter!(ax2, otu2.ax2.scatterx, otu2.ax2.scattery, color=:white, strokecolor=colors[2], markersize=4, strokewidth=0.4)
    scatter!(ax2, otu3.ax2.scatterx, otu3.ax2.scattery, color=:white, strokecolor=colors[3], markersize=4, strokewidth=0.4)
    scatter!(ax2, otu5.ax2.scatterx, otu5.ax2.scattery, color=:white, strokecolor=colors[4], markersize=4, strokewidth=0.4)
    lines!(ax2, otu5.ax2.linex, otu5.ax2.liney,
            linewidth=1, linestyle=:dash, color=:black)

    # Plot ax3
    scatter!(ax3, otu.ax3.scatterx, otu.ax3.scattery, color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4)
    scatter!(ax3, otu2.ax3.scatterx, otu2.ax3.scattery, color=:white, strokecolor=colors[2], markersize=4, strokewidth=0.4)
    scatter!(ax3, otu3.ax3.scatterx, otu3.ax3.scattery, color=:white, strokecolor=colors[3], markersize=4, strokewidth=0.4)
    scatter!(ax3, otu5.ax3.scatterx, otu5.ax3.scattery, color=:white, strokecolor=colors[4], markersize=4, strokewidth=0.4)
    lines!(ax3, otu5.ax3.linex, otu5.ax3.liney .^ 0.5,
            linewidth=1, linestyle=:dash, color=:black)

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
