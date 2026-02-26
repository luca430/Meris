#= Module to plot Taylor's law panels (flexible, embeddable like plot-SAD.jl) =#
#/ Start module
module TaylorPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings
using FileIO, ImageTransformations
using StatsBase, Random

using JLD2
using Colors, ColorTypes

#/ Modules
import Meris

const ICONDIR = Meris.FIGDIR .* "icons"
const MM_TO_PT = 72.0 / 25.4
const NATURE_SINGLE_WIDTH_PT = 89.0 * MM_TO_PT
const NATURE_DOUBLE_WIDTH_PT = 183.0 * MM_TO_PT
const NATURE_MAX_HEIGHT_PT = 170.0 * MM_TO_PT
const NATURE_AXIS_LABEL_PT = 7
const NATURE_TICK_PT = 6
const NATURE_TEXT_PT = 6
const NATURE_PANEL_LABEL_PT = 8
Random.seed!(1234) 

#################
### INTERNALS ###

"""Return a Vector of NamedTuples describing the default datasets."""
function _default_taylor_datasets(TLDIR)
    return [
        (;
            key=:linguistic,
            file=joinpath(TLDIR, "linguistic.jld2"),
            icon=joinpath(ICONDIR, "document.png"),
            icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.08, valign=0.92),
            occ_small=0.0,
            occ_big=0.9,
            take=5000,
            ref_shift=-0.5,
        ),
        (;
            key=:microbial,
            file=joinpath(TLDIR, "microbial.jld2"),
            icon=joinpath(ICONDIR, "bacteria.png"),
            icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.08,  valign=0.92),
            occ_small=0.0,
            occ_big=0.9,
            take=8000,
            ref_shift=-0.6,
        ),
        (;
            key=:social,
            file=joinpath(TLDIR, "social.jld2"),
            icon=joinpath(ICONDIR, "socio-economic.png"),
            icon_kw=(; width=Relative(0.77*0.25), height=Relative(0.25), halign=0.08,  valign=0.92),
            occ_small=0.0,
            occ_big=0.5,
            take=15000,
            ref_shift=-0.9,
        ),
        (;
            key=:biology,
            file=joinpath(TLDIR, "biology.jld2"),
            icon=joinpath(ICONDIR, "eco.png"),
            icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.08, valign=0.92),
            occ_small=0.0,
            occ_big=0.9,
            take=10000,
            ref_shift=-0.5,
        ),
    ]
end

"""Build a set of colors.

- If `palette` is provided, it is parsed and cycled if needed.
- Otherwise, it mimics plot-SAD.jl: base color from MakiePublication.COLORS, then HSL shades.
"""
function _make_colors(; palette, n::Int, color_num::Int=1, color_shades::Int=8)
    if !isnothing(palette)
        cols = parse.(Colorant, palette)
        return [cols[mod1(i, length(cols))] for i in 1:n]
    end

    base = MakiePublication.COLORS[1][color_num]
    base_hsl = HSL(base)
    shades = [HSL(base_hsl.h, base_hsl.s, l) for l in range(0.20, 0.80, length=max(color_shades, n))]
    return shades[1:n]
end

"""Load a JLD2 file and return the expected Taylor DataFrame under key "tldf"."""
function _load_tldf(path::AbstractString)
    return JLD2.load(path)["tldf"]
end

"""Compute centered log10 mean/var vectors (and drop nonpositive variance)."""
function _centered_logs(df; mean_col=:omeanfrequency, var_col=:ovarfrequency)
    dff = df[df[!, var_col] .> 0.0, :]
    m = log10.(dff[!, mean_col])
    s = log10.(dff[!, var_col])
    m .-= mean(m)
    s .-= mean(s)
    return dff, m, s
end

"""Add a small icon Axis into a given grid cell."""
function _add_icon!(parent_cell, icon_path;
    width=Relative(0.30), height=Relative(0.35),
    halign=0.07, valign=0.95
)
    axicon = Axis(
        parent_cell;
        width=width, height=height,
        halign=halign, valign=valign,
        tellwidth=false, tellheight=false,   # <-- this is the key
    )

    icon = FileIO.load(icon_path)
    icon_small = imresize(icon, (256, 256))
    image!(axicon, rotr90(icon_small))
    hidedecorations!(axicon)
    hidespines!(axicon)

    return axicon
end

"""Safe getter for NamedTuple fields."""
_ntget(nt::NamedTuple, field::Symbol, default) = hasproperty(nt, field) ? getproperty(nt, field) : default

#################
### FUNCTIONS ###

"""
    plot!(parent; kwargs...)

Create a Taylor's-law panel inside `parent` (a figure cell or GridLayout slot),
mirroring the embeddable style of `SADPlotter.plot!`.

Key knobs (similar spirit to plot-SAD.jl):
- `palette`: custom palette (hex strings) to control dataset colors.
- `color_num`, `color_shades`: fallback color generation when palette is `nothing`.
- `datasets`: override which datasets to plot (vector of NamedTuples like `_default_taylor_datasets`).
- `show_icons`: toggle the small icons.
- `big_limits`: limits for the big axis.
"""
function plot!(parent;
    TLDIR = Meris.DATADIR * "macro/taylor/",
    datasets = _default_taylor_datasets(TLDIR),
    palettes = [nothing, nothing, nothing, nothing],
    panel_start::Int = 1,
    font_scale::Float64 = 1.0,
    color_num = [1,2,3,4],
    color_shades = fill(8, 4),
    show_icons::Bool = true,
    big_limits = (-2, 2, -4, 4),
    small_limits = nothing,  # set to a 4-tuple, or keep auto
    panel_rowgap = 0,
    panel_colgap = 6,
    small_rowgap = 2,
    small_colgap = 3,
    markersize_small = 6,
    markersize_big = 6,
    strokewidth = 0.4,
)
    # Theme: keep the same MakiePublication look as SADPlotter
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    markers = [:circle, :rect, :diamond, :cross, :x, :utriangle, :dtriangle, :star4, :star6, :pentagon, :hexagon, :octagon]

    nsets = length(datasets)
    cols = Vector{Vector{Colorant}}(undef, nsets)

    # Layout inside the provided parent
    panel = GridLayout(parent)

    # Left: big axis spans two rows; Right: 2x2 grid of small axes
    ax_big = Axis(
        panel[1:2, 1],
        xlabel=L"\log_{10} \, \mu", xlabelsize=NATURE_AXIS_LABEL_PT * font_scale,
        ylabel=L"\log_{10} \, \sigma^2", ylabelsize=NATURE_AXIS_LABEL_PT * font_scale,
        xticklabelsize=NATURE_TICK_PT * font_scale, yticklabelsize=NATURE_TICK_PT * font_scale,
        aspect=AxisAspect(1),
        limits=big_limits,
    )

    right = GridLayout()
    panel[1:2, 2] = right

    rowsize!(panel, 1, Relative(0.5))
    rowsize!(panel, 2, Relative(0.5))
    colsize!(panel, 1, Relative(0.4))
    colsize!(panel, 2, Relative(0.55))
    rowgap!(panel, panel_rowgap)
    colgap!(panel, panel_colgap)
    rowgap!(right, small_rowgap)
    colgap!(right, small_colgap)

    # Create the 2x2 axes in the right grid
    axs_small = Axis[]
    grid_positions = [(2, 2), (2, 1), (1, 2), (1, 1)]
    for (i, (r, c)) in enumerate(grid_positions)
        i > nsets && break
        ax = Axis(
            right[r, c],
            xlabel=L"\log_{10} \, \mu", xlabelsize=NATURE_AXIS_LABEL_PT * font_scale,
            ylabel=L"\log_{10} \, \sigma^2", ylabelsize=NATURE_AXIS_LABEL_PT * font_scale,
            xticklabelsize=NATURE_TICK_PT * font_scale, yticklabelsize=NATURE_TICK_PT * font_scale,
        )
        if isnothing(small_limits)
            xlims!(ax, minimum(m_small_vec) * 1.5, maximum(m_small_vec) * 1.3)
            ylims!(ax, minimum(s_small_vec) * 1.5, maximum(s_small_vec) * 1.3)
        else
            xlims!(ax, small_limits[i][1], small_limits[i][2])
            ylims!(ax, small_limits[i][3], small_limits[i][4])
        end
        push!(axs_small, ax)
    end

    xtl = -10:0.1:10.0
    datasets = reverse(datasets)
    palettes = reverse(palettes)
    color_num = reverse(color_num)
    color_shades = reverse(color_shades)

    # Plot each dataset
    for (i, spec) in enumerate(datasets)
        i > length(axs_small) && break
        ax = axs_small[i]

        df = _load_tldf(spec.file)
        take = _ntget(spec, :take, :all)
        if take isa Integer
            total = size(df, 1)
            n = min(take, total)
            if n <= 0
                continue
            elseif n < total
                idx = randperm(total)[1:n]
                df = df[idx, :]
            end
        end

        classes = unique(df.class)
        cols[i] = _make_colors(;
            palette=palettes[i],
            n=length(classes),
            color_num=color_num[i],
            color_shades=color_shades[i]
        )
        class_order = reverse(classes)

        m_small_vec = Float64[]
        s_small_vec = Float64[]
        for (c,class) in enumerate(class_order)
            sdf = df[df.class .== class, :]
            df_small = spec.occ_small > 0 ? sdf[sdf.occupancy .> spec.occ_small, :] : sdf
            df_small, m_small, s_small = _centered_logs(df_small)
            isempty(m_small) && continue

            append!(m_small_vec, m_small)
            append!(s_small_vec, s_small)
            scatter!(
                ax, m_small, s_small;
                color=(:white, 1.0),
                strokecolor=cols[i][c],
                marker=markers[mod1(c, length(markers))],
                markersize=markersize_small,
                strokewidth=strokewidth,
            )
        end

        # same guide lines you used: y=2x and y=x, with the slight x-shift
        if !isempty(m_small_vec) && !isempty(s_small_vec)
            x0 = minimum(m_small_vec)
            y0 = minimum(s_small_vec)
            lines!(ax, xtl, 2 .* (xtl .- x0) .+ y0; linewidth=2, color=:black, linestyle=(:dash, :dense))
            lines!(ax, xtl, (xtl .- x0) .+ y0; linewidth=2, color=:grey, linestyle=(:dash, :dense))
        end

        # icon overlay
        if show_icons && isfile(spec.icon)
            ikw = _ntget(spec, :icon_kw, (;))  # default empty NamedTuple
            _add_icon!(right[grid_positions[i]...], spec.icon; ikw...)
        end

        # big axis: occupancy-filtered points from the same sampled dataset
        for (c,class) in enumerate(class_order)
            sdf = df[df.class .== class, :]
            df_big = spec.occ_big > 0 ? sdf[sdf.occupancy .> spec.occ_big, :] : sdf
            df_big, m_big, s_big = _centered_logs(df_big)
            isempty(m_big) && continue

            scatter!(
                ax_big, m_big, s_big;
                color=(:white, 1.0),
                strokecolor=cols[i][c],
                marker=markers[mod1(c, length(markers))],
                markersize=markersize_big,
                strokewidth=strokewidth,
            )
        end
    end

    # Reference lines in all axes
    lines!(ax_big, xtl, 2 .* xtl; linewidth=2, color=:black, linestyle=(:dash, :dense))

    # Panel letters outside axes, following Makie layout-label style.
    letters = [Char(Int('a') + mod(panel_start - 1 + i, 26)) for i in 0:4]
    Label(panel[1, 1, TopRight()], string(letters[1]);
        fontsize=NATURE_PANEL_LABEL_PT * font_scale, font=:bold, color=:black,
        halign=:right, valign=:bottom, padding=(0, 6, 6, 0)
    )
    label_positions = [(1, 1), (1, 2), (2, 1), (2, 2)]  # visual order: TL, TR, BL, BR
    for (i, (r, c)) in enumerate(label_positions)
        i > length(axs_small) && break
        Label(right[r, c, TopLeft()], string(letters[i + 1]);
            fontsize=NATURE_PANEL_LABEL_PT * font_scale, font=:bold, color=:black,
            halign=:right, valign=:bottom, padding=(0, 6, 6, 0)
        )
    end

    return (;
        panel,
        ax_big,
        axs_small,
        colors=cols,
    )
end

"""Backward-compatible wrapper that creates a standalone Figure."""
function plot_taylor(; TLDIR=Meris.DATADIR * "macro/taylor/", savefig=false, figname="tl.png", kwargs...)
    fig = Figure(; size=(NATURE_DOUBLE_WIDTH_PT, 0.65 * NATURE_MAX_HEIGHT_PT), figure_padding=(4, 4, 4, 4))

    plot!(fig[1, 1]; TLDIR=TLDIR, kwargs...)

    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module
