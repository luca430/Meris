#= Module to plot Zipf, CCDF and SAD considering a filter to remove low counts =#
#/ Start module
module SADPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta
using StatsBase, JLD2
using Colors, ColorTypes
using LsqFit

#/ Modules
import Meris
import Meris.MDistributions as MDist

#################
### FUNCTIONS ###
function plot!(parent;
        color_num=1,
        color_shades=4,
        palette=nothing,
        ZIPFDIR=Meris.DATADIR * "macro/sad/",
    )

    # Theme (you may want to set this ONCE outside when doing 4 panels)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    markers = [:circle, :rect, :diamond, :cross, :x, :utriangle, :dtriangle, :star4, :star6, :pentagon, :hexagon, :octagon]

    base = MakiePublication.COLORS[1][color_num]
    base_hsl = HSL(base)
    colors = [HSL(base_hsl.h, base_hsl.s, l) for l in range(0.15, 0.75, length=color_shades)]
    colors = isnothing(palette) ? colors : parse.(Colorant, palette)

    # ---- THIS is the key change: build a panel layout inside `parent`
    panel = GridLayout(parent)

    # Top + bottom sublayouts
    top = GridLayout()
    panel[1, 1] = top

    bottom = GridLayout()
    panel[2, 1] = bottom

    ax1 = Axis(top[1,1],
        xticklabelsize=6,
        xaxisposition = :top,
        xlabel="", ylabel=L"\text{CAD exponent } \gamma",
        ylabelsize=9,
        yminorticksvisible = false,
        xminorticksvisible = false,
        xticksvisible = false,
        ygridvisible=true,
        ygridwidth = 0.3
    )

    ax2 = Axis(
        bottom[1,1],
        xlabel=L"\text{rel. abundance } \nu",
        ylabel=L"p(\nu)",
        xlabelsize=9,
        ylabelsize=9,
        xscale=log10, yscale=log10,
    )

    ax3 = Axis(
        bottom[1,2],
        xlabel = L"\text{sample size } N",
        ylabel = L"\eta = \gamma - 1",
        xlabelsize = 9,
        ylabelsize = 9,
        xscale = log10, #yscale=log10,
        yticksmirrored = true
    )

    # sizes INSIDE the panel (not the global figure)
    rowsize!(panel, 1, Relative(0.4))
    rowsize!(panel, 2, Relative(0.6))
    rowgap!(panel, 15)

    colsize!(top, 1, Relative(0.91))

    colsize!(bottom, 1, Relative(0.45))
    colsize!(bottom, 2, Relative(0.45))

    # ---- Data + plotting (unchanged)
    out = JLD2.load(ZIPFDIR)["out"]

    labels = collect(keys(out.pl))
    pl_vals = collect(values(out.pl))

    min_vals = Float64[]
    max_vals = Float64[]
    for (i, p) in enumerate(pl_vals)
        push!(min_vals, minimum(p.α) + 1)
        push!(max_vals, maximum(p.α) + 1)
        boxplot!(ax1, fill(i, length(p.α)), p.α .+ 1,
            color=(colors[i], 0.8), markersize=4, whiskerlinewidth=0.8,
            medianlinewidth=0.6, width=0.4)
    end

    ymin, ymax = minimum(min_vals), maximum(max_vals)
    ylims!(ax1, ymin, ymax)
    ypos = ymax - 0.1*(ymax - ymin)

    for i in 1:length(labels)
        scatter!(ax1, [i], [ypos];
            marker = markers[i],
            color = :white,
            strokecolor = colors[i],
            markersize = 8,
            strokewidth = 0.5
        )
    end
    ax1.xticks = (1:length(labels), labels)

    cad_vals = collect(values(out.cad))
    x_min = Float64[]
    x_max = Float64[]
    for (i, (v, p)) in enumerate(zip(cad_vals, pl_vals))
        push!(x_min, minimum(v.x))
        push!(x_max, maximum(v.x))
        α = p.α_eff
        ε = p.ε_eff
        scatter!(ax2, 10 .^ (v.x), (v.y ./ (α * ε ^ α .* log(10))) .^ (1 / α),
            marker=markers[i], color=:white, strokecolor=colors[i],
            markersize=5, strokewidth=0.4)
    end
    xrange = minimum(x_min)*1.05:1e-2:maximum(x_max)/1.1
    lines!(ax2, 10 .^ xrange, 10 .^ (-xrange), color=:black, linestyle=:dash, linewidth=1)

    heaps_vals = collect(values(out.heaps))
    for (i, (v,p)) in enumerate(zip(heaps_vals, pl_vals))
        α = p.α_eff
        z, t = log10.(v.N), log10.(v.V)
        g = diff(t) ./ diff(z)
        c = (g[end] - 1) / (α - 1)
        etas = (g .- 1) ./ c .+ 1
        idx = round.(Int, 10 .^ range(0, log10(length(etas)), length = 60))
        scatter!(ax3, v.N[idx], etas[idx],
            marker=markers[i], color=:white, strokecolor=colors[i],
            markersize=5, strokewidth=0.4)
            hlines!(ax3, α, color=colors[i], linestyle=:dot, linewidth=0.5, xmin=(log10(v.N[idx][end])-1) / 6)
    end
    hlines!(ax3, 1, color=:black, linestyle=:dash, linewidth=0.5)

    return (ax1=ax1, ax2=ax2, ax3=ax3, panel=panel)
end

function plot(; savefig=false, figname="sad.png", kwargs...)
    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(1.75*width, 1.5*height), figure_padding=(2,4,2,14))
    plot!(fig[1,1]; kwargs...)
    (savefig && !isnothing(figname)) && CairoMakie.save(figname, fig, pt_per_unit=1)
    return fig
end


### HELPER ###
function plot4(;
        ZIPFDIR=Meris.DATADIR * "macro/sad/",
        savefig=false,
        figname="full_sad.png",
        pt_per_unit=1
    )

    width = .95 * 246
    height = 3*width / 4.67

    bigfig = Figure(
        size = (2 * 1.75*width, 2 * 1.5*height),
        figure_padding = (8, 8, 8, 8)
    )

    # 2×2 panels
    plot!(bigfig[1,1]; ZIPFDIR=ZIPFDIR, color_num=1)
    plot!(bigfig[1,2]; ZIPFDIR=ZIPFDIR, color_num=2)
    plot!(bigfig[2,1]; ZIPFDIR=ZIPFDIR, color_num=3)
    plot!(bigfig[2,2]; ZIPFDIR=ZIPFDIR, color_num=4)

    rowgap!(bigfig.layout, 20)
    colgap!(bigfig.layout, 20)

    savefig && CairoMakie.save(figname, bigfig, pt_per_unit=pt_per_unit)
    return bigfig
end

end # module SADPlotter
#/ End module
