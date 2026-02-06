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
        ax1limits=(nothing, nothing, nothing, nothing),
        ax2limits=(nothing, nothing, nothing, nothing),
        ax3limits=(nothing, nothing, nothing, nothing)
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
        ygridwidth = 0.3,
        limits=ax1limits
    )

    ax2 = Axis(
        bottom[1,1],
        xlabel=L"\text{rel. abundance } \nu",
        ylabel=L"p(\nu)",
        xlabelsize=9,
        ylabelsize=9,
        xscale=log10, yscale=log10,
        limits=ax2limits
    )

    ax3 = Axis(
        bottom[1,2],
        xlabel = L"\text{sample size } N",
        ylabel = L"\text{sample size } N",
        xlabelsize = 9,
        ylabelsize = 9,
        xscale = log10, yscale=log10,
        limits=ax3limits
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

    ### AX1 ###
    labels = sort(collect(keys(out.pl)))
    pl_vals = [out.pl[k] for k in labels]

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

    ### AX2 ###
    cad_vals = [out.cad[k] for k in labels]
    ple_vals = [out.ple[k] for k in labels]
    x_min = Float64[]
    x_max = Float64[]
    for (i, (v, p, pe)) in enumerate(zip(cad_vals, pl_vals, ple_vals))
        push!(x_min, minimum(v.x))
        push!(x_max, maximum(v.x))
        scatter!(ax2, 10 .^ (v.x), (v.y ./ (pe.α[1] * mean(p.ε) ^ pe.α[1] .* log(10))) .^ (1 / pe.α[1]),
            marker=markers[i], color=:white, strokecolor=colors[i],
            markersize=5, strokewidth=0.4)
    end
    xrange = minimum(x_min)*1.05:1e-2:maximum(x_max)/1.1
    lines!(ax2, 10 .^ xrange, 10 .^ (-xrange), color=:black, linestyle=:dash, linewidth=1.5)

    ### AX3 ###
    heaps_vals = [out.heaps[k] for k in labels]
    α_vals = Float64[]
    xmax = []
    for (i, (v,p)) in enumerate(zip(heaps_vals, pl_vals))
        push!(α_vals, mean(p.α))
        push!(xmax, maximum(v.N))
        idx = round.(Int, 10 .^ range(0, log10(length(v.V)), length = 60))
        scatter!(ax3, v.N[idx], v.V[idx],
            marker=markers[i], color=:white, strokecolor=colors[i],
            markersize=5, strokewidth=0.4)
    end

    parss = []
    for (i, v) in enumerate(heaps_vals)
        a = minimum([α_vals[i], 1])
        fit, pars, which = Meris.HeapsModel.fit_regimes(v.N, v.V; a=a)
        push!(parss, pars)
    end

    text!(ax3, 5e1, 5e1 * 2, text=L"V \sim N", color=colors[1], rotation = π/4, fontsize=10)
    xp = maximum(xmax) / 500
    yp = Meris.HeapsModel.predict_regimes(xp, parss[argmax(α_vals)])
    text!(ax3, xp , yp * 1.5, text=L"V \sim N^{\eta}", color=colors[1], rotation = atan(minimum([maximum(α_vals),1])), fontsize=10)

    # --- INSET in ax3: zoom on the tail (end of curves) ---
    # place a small axis on top of ax3 (relative to ax3 cell)
    inset = Axis(
        bottom[1,2],
        width  = Relative(0.43),
        height = Relative(0.43),
        halign = 0.95,   # right
        valign = 0.12,   # bottom
        xscale = log10,
        yscale = log10,
        xticklabelsize = 5,
        yticklabelsize = 5,
        xlabelsize = 7,
        ylabelsize = 7,
        xgridvisible = false,
        ygridvisible = false,
        yminorticksvisible = false,
        xminorticksvisible = false,
        limits=ax3limits
    )

    # choose what "end" means (last decade by default)
    x_hi = maximum(xmax)
    x_lo = 10            # last decade; change to /30, /100 if you want

    # compute corresponding y-lims from the fitted curves (robust and consistent)
    yy = Float64[]
    for (i, v) in enumerate(heaps_vals)
        xr = 10 .^ range(log10(x_lo), log10(x_hi); length=200)
        yhat = Meris.HeapsModel.predict_regimes(xr, parss[i])
        append!(yy, yhat)
    end

    # replot tail points + tail fit lines into the inset
    for (i, v) in enumerate(heaps_vals)
        # tail indices
        idx_tail = findall(n -> n ≥ x_lo, v.N)
        xr = 10 .^ range(log10(x_lo), log10(x_hi); length=200)
        yhat = Meris.HeapsModel.predict_regimes(xr, parss[i])
        lines!(inset, xr, yhat, color=colors[i], linestyle=:dash, linewidth=0.7)
    end

    # optional: hide labels (usually nicer for insets)
    inset.xlabel = ""
    inset.ylabel = ""
    inset.xticklabelsvisible = false
    inset.yticklabelsvisible = false

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
