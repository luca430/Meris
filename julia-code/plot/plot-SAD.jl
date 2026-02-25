#= Module to plot Zipf, CCDF and SAD considering a filter to remove low counts =#
#/ Start module
module SADPlotter

using CairoMakie
using Makie: Point2f, apply_transform, transformationmatrix
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta
using StatsBase, JLD2
using Colors, ColorTypes
using LsqFit
using FileIO, ImageTransformations

#/ Modules
import Meris
import Meris.MDistributions as MDist

const ICONDIR = Meris.FIGDIR .* "icons"

#################
### FUNCTIONS ###
function plot!(
    parent;
    color_num=1,
    color_shades=4,
    palette=nothing,
    ZIPFDIR=Meris.DATADIR * "macro/sad/",
    ax1limits=(nothing, nothing, nothing, nothing),
    ax2limits=(nothing, nothing, nothing, nothing),
    ax3limits=(nothing, nothing, nothing, nothing),
    reverse_panel=false,
    icon_name=nothing,
    icon_kw=(; width=Relative(0.25), height=Relative(0.3), halign=0.05, valign=0.95),
    ax2_text_offset=(1.0, 4.0),
    ax3_text_color=nothing,
    ax3_text_offset=(1.0, 1.5)
    )

    # Theme
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    markers = [:circle, :rect, :diamond, :cross, :x, :utriangle, :dtriangle, :star4, :star6, :pentagon, :hexagon, :octagon]

    base = MakiePublication.COLORS[1][color_num]
    base_hsl = HSL(base)
    colors = [HSL(base_hsl.h, base_hsl.s, l) for l in range(0.15, 0.75, length=color_shades)]
    colors = isnothing(palette) ? colors : parse.(Colorant, palette)

    # Build a panel layout inside `parent`
    panel = GridLayout(parent)

    # Top + bottom sublayouts
    top_row = reverse_panel ? 2 : 1
    bottom_row = reverse_panel ? 1 : 2

    top = GridLayout()
    panel[top_row, 1] = top

    bottom = GridLayout()
    panel[bottom_row, 1] = bottom

    ax1_cell = top[1,1]
    ax2_cell = bottom[1,1]
    ax3_cell = bottom[1,2]

    ax1 = Axis(ax1_cell,
        xaxisposition = reverse_panel ? :bottom : :top,
        xlabel="", ylabel=L"\text{CAD } \gamma",
        ylabelsize=12,
        xticklabelalign = (:right, :center),
        xticklabelrotation = reverse_panel ? π/6 : -π/6,
        xticklabelpad = 5,
        xticklabelsize=9,
        yticklabelsize=10,
        yminorticksvisible = false,
        xminorticksvisible = false,
        xticksvisible = false,
        ygridvisible=true,
        ygridwidth = 0.3,
        limits=ax1limits
    )

    ax2 = Axis(
        ax2_cell,
        xlabel=L"\text{rel. abundance } \nu",
        ylabel=L"p(\nu)",
        xlabelsize=12,
        ylabelsize=12,
        xticklabelsize=10,
        yticklabelsize=10,
        xscale=log10, yscale=log10,
        limits=ax2limits
    )

    ax3 = Axis(
        ax3_cell,
        xlabel = L"\text{sample size } N",
        ylabel = L"\text{vocabulary size } V",
        xlabelsize = 12,
        ylabelsize = 12,
        xticklabelsize=10,
        yticklabelsize=10,
        xscale = log10, yscale=log10,
        limits=ax3limits
    )

    # sizes INSIDE the panel (not the global figure)
    if reverse_panel
        rowsize!(panel, 1, Relative(0.6))
        rowsize!(panel, 2, Relative(0.4))
    else
        rowsize!(panel, 1, Relative(0.4))
        rowsize!(panel, 2, Relative(0.6))
    end
    rowgap!(panel, 15)

    colsize!(top, 1, Relative(0.91))

    colsize!(bottom, 1, Relative(0.45))
    colsize!(bottom, 2, Relative(0.45))

    # Data + plotting
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

    data_ymin = minimum(min_vals)
    data_ymax = maximum(max_vals)

    user_ymin = ax1limits[3]
    user_ymax = ax1limits[4]
    ymin = isnothing(user_ymin) ? data_ymin * 0.8 : user_ymin
    ymax = isnothing(user_ymax) ? data_ymax * 1.2 : user_ymax
    if ymax <= ymin
        ymax = ymin + max(abs(ymin) * 0.1, 1e-3)
    end
    ylims!(ax1, ymin, ymax)

    ax1_lims = ax1.limits[]
    ax1_ymin = (ax1_lims isa Tuple && length(ax1_lims) >= 4 && !isnothing(ax1_lims[3])) ? ax1_lims[3] : ymin
    ax1_ymax = (ax1_lims isa Tuple && length(ax1_lims) >= 4 && !isnothing(ax1_lims[4])) ? ax1_lims[4] : ymax
    span = ax1_ymax - ax1_ymin
    ypos = reverse_panel ? (ax1_ymin + 0.08 * span) : (ax1_ymax - 0.08 * span)

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

    ax2_lims = ax2.limits[]
    ax2_xmin = (ax2_lims isa Tuple && length(ax2_lims) >= 4 && !isnothing(ax2_lims[1])) ? ax2_lims[1] : minimum(x_min)
    ax2_xmax = (ax2_lims isa Tuple && length(ax2_lims) >= 4 && !isnothing(ax2_lims[2])) ? ax2_lims[2] : maximum(x_max)
    xrange = (ax2_xmin-1.1):1e-2:-0.5
    lines!(ax2, 10 .^ xrange, 10 .^ (-xrange), color=:black, linestyle=:dash, linewidth=1.5, label=L"y \sim x^{-1}")
    axislegend(ax2,
        position=:rt,
        labelsize=12,
        markersize=2,
        patchsize=(12, 5)
    )

    if !isnothing(icon_name)
        icon_path = joinpath(ICONDIR, icon_name)
        add_icon!(ax2_cell, icon_path; icon_kw...)
    end

    ### AX3 ###
    heaps_vals = [out.heaps[k] for k in labels]
    α_vals = Float64[]
    xmax = []
    ymax = []
    for (i, (v,p)) in enumerate(zip(heaps_vals, pl_vals))
        push!(α_vals, mean(p.α))
        push!(xmax, maximum(v.N))
        push!(ymax, maximum(v.V))
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

    ax3_text_color = isnothing(ax3_text_color) ? colors[1] : ax3_text_color
    text!(ax3, 5e1, 5e1 * 2, text=L"V \sim N", color=ax3_text_color, rotation = π/4, fontsize=10)

    # pick two nearby x points on the curve (in DATA coordinates): this is because we need to express data slope in terms of axis dimensions
    xp = maximum(xmax) / 500
    x1 = xp
    x2 = xp * 1.15
    
    p  = parss[argmax(ymax)]
    y1 = Meris.HeapsModel.predict_regimes(x1, p)
    y2 = Meris.HeapsModel.predict_regimes(x2, p)
    
    M = Makie.transformationmatrix(ax3.scene)[]
    v1 = Vec4f(Float32(x1), Float32(y1), 0f0, 1f0)
    v2 = Vec4f(Float32(x2), Float32(y2), 0f0, 1f0)
    
    p1 = M * v1
    p2 = M * v2
    
    θ = atan(p2[2] - p1[2], p2[1] - p1[1])  # atan(dy, dx)
    
    text!(ax3, x1 * ax3_text_offset[1], y1 * ax3_text_offset[2];
        text=L"V \sim N^{\eta}",
        rotation=θ,
        color=ax3_text_color,
        fontsize=10
    )

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
function add_icon!(cell, icon_path;
    width=Relative(0.20), height=Relative(0.20),
    halign=0.90, valign=0.90,
    rotate=true, pixels=256
)
    axicon = Axis(cell;
        width=width, height=height,
        halign=halign, valign=valign,
        tellwidth=false, tellheight=false,
    )

    img = FileIO.load(icon_path)
    img = imresize(img, (pixels, pixels))
    rotate && (img = rotr90(img))

    image!(axicon, img)
    hidedecorations!(axicon)
    hidespines!(axicon)
    return axicon
end

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
