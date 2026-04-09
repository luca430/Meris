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
const MM_TO_PT = 72.0 / 25.4
const NATURE_SINGLE_WIDTH_PT = 89.0 * MM_TO_PT
const NATURE_DOUBLE_WIDTH_PT = 183.0 * MM_TO_PT
const NATURE_MAX_HEIGHT_PT = 170.0 * MM_TO_PT
const NATURE_AXIS_LABEL_PT = 7
const NATURE_TICK_PT = 6
const NATURE_TEXT_PT = 6
const NATURE_PANEL_LABEL_PT = 8

#################
### FUNCTIONS ###
function plot!(parent;
        color_num=1,
        panel_id=1,
        font_scale=1.0,
        color_shades=4,
        palette=nothing,
        DIR=Meris.DATADIR * "macro/sad/",
        CADDIR=nothing,
        ax1limits=(nothing, nothing, nothing, nothing),
        ax2limits=(nothing, nothing, nothing, nothing),
        ax3limits=(nothing, nothing, nothing, nothing),
        reverse_panel=false,
        icon_name=nothing,
        icon_kw=(; width=Relative(0.25), height=Relative(0.3), halign=0.05, valign=0.95),
        ax2_text_offset=(1.0, 0.2),
        ax3_text_offset=(1.0, 1.5),
        inner_rowgap=3
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
        ylabelsize=NATURE_AXIS_LABEL_PT * font_scale,
        xticklabelalign = (:right, :center),
        xticklabelrotation = reverse_panel ? π/6 : -π/6,
        xticklabelpad = 5,
        xticklabelsize=NATURE_TICK_PT * font_scale,
        yticklabelsize=NATURE_TICK_PT * font_scale,
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
        xlabelsize=NATURE_AXIS_LABEL_PT * font_scale,
        ylabelsize=NATURE_AXIS_LABEL_PT * font_scale,
        xticklabelsize=NATURE_TICK_PT * font_scale,
        yticklabelsize=NATURE_TICK_PT * font_scale,
        xscale=log10, yscale=log10,
        limits=ax2limits
    )

    ax3 = Axis(
        ax3_cell,
        xlabel = L"\text{sample size } N",
        ylabel = L"\text{vocabulary size } V",
        xlabelsize = NATURE_AXIS_LABEL_PT * font_scale,
        ylabelsize = NATURE_AXIS_LABEL_PT * font_scale,
        xticklabelsize=NATURE_TICK_PT * font_scale,
        yticklabelsize=NATURE_TICK_PT * font_scale,
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
    rowgap!(panel, inner_rowgap)

    colsize!(top, 1, Relative(0.91))

    colsize!(bottom, 1, Relative(0.45))
    colsize!(bottom, 2, Relative(0.45))

    # Data + plotting
    target = isnothing(CADDIR) ? DIR : CADDIR
    data_root = isdir(target) ? target : dirname(target)
    files = if isdir(target)
        sort(filter(f -> endswith(lowercase(f), ".jld2"), readdir(target)))
    elseif isfile(target)
        [basename(target)]
    else
        String[]
    end
    isempty(files) && error("No .jld2 files found in $(target)")

    series = NamedTuple[]
    for file in files
        out = JLD2.load(joinpath(data_root, file))["out"]
        labels = sort(collect(keys(out.pl)))
        for label in labels
            haskey(out.cad, label) || continue
            haskey(out.heaps, label) || continue
            push!(series, (;
                label = string(label),
                file = file,
                pl = out.pl[label],
                cad = out.cad[label],
                heaps = out.heaps[label]
            ))
        end
    end
    isempty(series) && error("No pl/cad/heaps series found in $(target)")

    raw_labels = [s.label for s in series]
    xticklabels = if length(unique(raw_labels)) == length(raw_labels)
        raw_labels
    else
        ["$(s.label)\n($(splitext(s.file)[1]))" for s in series]
    end

    min_vals = Float64[]
    max_vals = Float64[]
    x_min = Float64[]
    x_max = Float64[]
    α_vals = Float64[]
    xmax = Float64[]
    ymax = Float64[]
    parss = Any[]

    for (i, s) in enumerate(series)
        c = colors[mod1(i, length(colors))]
        m = markers[mod1(i, length(markers))]

        p = s.pl
        push!(min_vals, minimum(p.α) + 1)
        push!(max_vals, maximum(p.α) + 1)
        boxplot!(ax1, fill(i, length(p.α)), p.α .+ 1,
            color=(c, 0.8), markersize=4, whiskerlinewidth=0.8,
            medianlinewidth=0.6, width=0.4)

        v = s.cad
        push!(x_min, minimum(v.x))
        push!(x_max, maximum(v.x))
        scatter!(ax2, 10 .^ (v.x), (v.y ./ (mean(p.α) * mean(p.ε) ^ mean(p.α) .* log(10))) .^ (1 / mean(p.α)),
            marker=m, color=:white, strokecolor=c,
            markersize=5, strokewidth=0.4)

        h = s.heaps
        push!(α_vals, mean(p.α))
        push!(xmax, maximum(h.N))
        push!(ymax, maximum(h.V))
        idx = round.(Int, 10 .^ range(0, log10(length(h.V)), length = 60))
        scatter!(ax3, h.N[idx], h.V[idx],
            marker=m, color=:white, strokecolor=c,
            markersize=5, strokewidth=0.4)

        a = minimum([α_vals[end], 1])
        fit, pars, which = Meris.HeapsModel.fit_regimes(h.N, h.V; a=a)
        push!(parss, pars)
    end

    ### AX1 ###
    data_ymin = minimum(min_vals)
    data_ymax = maximum(max_vals)
    user_ymin = ax1limits[3]
    user_ymax = ax1limits[4]
    ymin = isnothing(user_ymin) ? data_ymin * 0.8 : user_ymin
    ymax_ax1 = isnothing(user_ymax) ? data_ymax * 1.2 : user_ymax
    if ymax_ax1 <= ymin
        ymax_ax1 = ymin + max(abs(ymin) * 0.1, 1e-3)
    end
    ylims!(ax1, ymin, ymax_ax1)

    ax1_lims = ax1.limits[]
    ax1_ymin = (ax1_lims isa Tuple && length(ax1_lims) >= 4 && !isnothing(ax1_lims[3])) ? ax1_lims[3] : ymin
    ax1_ymax = (ax1_lims isa Tuple && length(ax1_lims) >= 4 && !isnothing(ax1_lims[4])) ? ax1_lims[4] : ymax_ax1
    span = ax1_ymax - ax1_ymin
    ypos = reverse_panel ? (ax1_ymin + 0.08 * span) : (ax1_ymax - 0.08 * span)

    for i in 1:length(series)
        c = colors[mod1(i, length(colors))]
        m = markers[mod1(i, length(markers))]
        scatter!(ax1, [i], [ypos];
            marker = m,
            color = :white,
            strokecolor = c,
            markersize = 8,
            strokewidth = 0.5
        )
    end
    ax1.xticks = (1:length(series), xticklabels)

    angle_in_screen(ax, x1, y1, x2, y2) = begin
        tf = ax.scene.transformation.transform_func[]
        q1 = apply_transform(tf, Point2f(Float32(x1), Float32(y1)))
        q2 = apply_transform(tf, Point2f(Float32(x2), Float32(y2)))
        M = Makie.transformationmatrix(ax.scene)[]
        v1 = [Float64(q1[1]), Float64(q1[2]), 0.0, 1.0]
        v2 = [Float64(q2[1]), Float64(q2[2]), 0.0, 1.0]
        p1 = M * v1
        p2 = M * v2
        atan(p2[2] - p1[2], p2[1] - p1[1])  # atan(dy, dx)
    end

    ### AX2 ###
    ax2_lims = ax2.limits[]
    ax2_xmin = (ax2_lims isa Tuple && length(ax2_lims) >= 4 && !isnothing(ax2_lims[1])) ? ax2_lims[1] : minimum(x_min)
    ax2_xmax = (ax2_lims isa Tuple && length(ax2_lims) >= 4 && !isnothing(ax2_lims[2])) ? ax2_lims[2] : maximum(x_max)
    xrange = log10(ax2_xmin):1e-2:log10(ax2_xmax)
    lines!(ax2, 10 .^ xrange, 10 .^ (-xrange), color=:black, linestyle=:dash, linewidth=1.5)
    x2_ref = sqrt(ax2_xmin * ax2_xmax)
    x2_ref2 = x2_ref * 1.15
    y2_ref = x2_ref^-1
    y2_ref2 = x2_ref2^-1
    θ2 = angle_in_screen(ax2, x2_ref, y2_ref, x2_ref2, y2_ref2)
    text!(ax2, x2_ref * ax2_text_offset[1], y2_ref * ax2_text_offset[2];
        text=L"\sim \nu^{-1}",
        rotation=θ2,
        color=:black,
        fontsize=10
    )

    if !isnothing(icon_name)
        icon_path = joinpath(ICONDIR, icon_name)
        add_icon!(ax2_cell, icon_path; icon_kw...)
    end

    ### AX3 ###
    ax3_lims = ax3.limits[]
    ax3_xmin = (ax3_lims isa Tuple && length(ax3_lims) >= 4 && !isnothing(ax3_lims[1])) ? ax3_lims[1] : minimum([minimum(s.heaps.N) for s in series])
    ax3_xmax = (ax3_lims isa Tuple && length(ax3_lims) >= 4 && !isnothing(ax3_lims[2])) ? ax3_lims[2] : maximum(xmax)
    x3_ref = 50
    x3_ref2 = x3_ref * 1.15
    y3_ref = x3_ref
    y3_ref2 = x3_ref2
    θ3_id = angle_in_screen(ax3, x3_ref, y3_ref, x3_ref2, y3_ref2)
    text!(ax3, x3_ref * ax3_text_offset[1], y3_ref * ax3_text_offset[2];
        text=L"V \sim N",
        color=:black,
        rotation=θ3_id,
        fontsize=9
    )

    # pick two nearby x points on the curve (in DATA coordinates): this is because we need to express data slope in terms of axis dimensions
    xp = ax3_xmax / 400
    x1 = xp
    x2 = xp * 1.15

    p = parss[argmax(ymax)]
    y1 = Meris.HeapsModel.predict_regimes(x1, p)
    y2 = Meris.HeapsModel.predict_regimes(x2, p)
    θ = angle_in_screen(ax3, x1, y1, x2, y2)

    text!(ax3, x1 * ax3_text_offset[1], y1 * ax3_text_offset[2];
        text=L"V \sim N^{\eta}",
        rotation=θ,
        color=:black,
        fontsize=9
    )

    base_letter_idx = 3 * (panel_id - 1)
    letters = [
        Char(Int('A') + mod(base_letter_idx + 0, 26)),
        Char(Int('A') + mod(base_letter_idx + 1, 26)),
        Char(Int('A') + mod(base_letter_idx + 2, 26))
    ]

    # Panel letters outside axes, following Makie layout-label style.
    Label(top[1,1,TopRight()], string(letters[1]);
        fontsize=NATURE_PANEL_LABEL_PT * font_scale, font=:bold, color=:black,
        halign=:right, valign=:bottom, padding=(0, 6, 6, 0)
    )
    Label(bottom[1,1,TopRight()], string(letters[2]);
        fontsize=NATURE_PANEL_LABEL_PT * font_scale, font=:bold, color=:black,
        halign=:right, valign=:bottom, padding=(0, 6, 6, 0)
    )
    Label(bottom[1,2,TopRight()], string(letters[3]);
        fontsize=NATURE_PANEL_LABEL_PT * font_scale, font=:bold, color=:black,
        halign=:right, valign=:bottom, padding=(0, 6, 6, 0)
    )

    return (ax1=ax1, ax2=ax2, ax3=ax3, panel=panel)
end

function plot(; savefig=false, figname="sad.png", kwargs...)
    fig = Figure(;
        size=(NATURE_SINGLE_WIDTH_PT, 0.55 * NATURE_SINGLE_WIDTH_PT),
        figure_padding=(4, 4, 4, 4)
    )
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
        CADDIR=Meris.DATADIR * "macro/sad/",
        savefig=false,
        figname="full_sad.png",
        pt_per_unit=1
    )
    bigfig = Figure(
        size = (NATURE_DOUBLE_WIDTH_PT, NATURE_MAX_HEIGHT_PT),
        figure_padding = (4, 4, 4, 4)
    )

    # 2×2 panels
    plot!(bigfig[1,1]; CADDIR=CADDIR, color_num=1, panel_id=1)
    plot!(bigfig[1,2]; CADDIR=CADDIR, color_num=2, panel_id=2)
    plot!(bigfig[2,1]; CADDIR=CADDIR, color_num=3, panel_id=3)
    plot!(bigfig[2,2]; CADDIR=CADDIR, color_num=4, panel_id=4)

    rowgap!(bigfig.layout, 20)
    colgap!(bigfig.layout, 20)

    savefig && CairoMakie.save(figname, bigfig, pt_per_unit=pt_per_unit)
    return bigfig
end

end # module SADPlotter
#/ End module
