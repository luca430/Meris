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

const ICONDIR = normpath(joinpath(@__DIR__, "..", "icons"))
Random.seed!(1234) 

#################
### INTERNALS ###

"""Return a Vector of NamedTuples describing the default datasets."""
function _default_taylor_datasets(TLDIR)
    return [
        (;
            key=:linguistic,
            file=joinpath(TLDIR, "linguistic.jld2"),
            icon=joinpath(ICONDIR, "documents.png"),
            occ_small=0.0,
            occ_big=0.9,
            take=400,
            # small panels: shift the reference lines a bit (kept from your original)
            ref_shift=-0.5,
        ),
        (;
            key=:microbial,
            file=joinpath(TLDIR, "microbial.jld2"),
            icon=joinpath(ICONDIR, "bacteria.png"),
            occ_small=0.0,   # show all in small panel
            occ_big=0.9,
            take=200,
            ref_shift=-0.6,
        ),
        (;
            key=:social,
            file=joinpath(TLDIR, "social.jld2"),
            icon=joinpath(ICONDIR, "lego.png"),
            occ_small=0.0,
            occ_big=0.5,
            take=400,
            ref_shift=-0.9,
        ),
        (;
            key=:biology,
            file=joinpath(TLDIR, "biology.jld2"),
            icon=joinpath(ICONDIR, "gene.png"),
            occ_small=0.0,
            occ_big=0.99999,
            take=400,
            ref_shift=-1.3,
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
        return [cols[mod1(i, length(cols))] for i in 1:length(cols)]
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
function _add_icon!(parent_cell, icon_path; width=Relative(0.22), height=Relative(0.22), halign=0.072, valign=0.95)
    axicon = Axis(parent_cell; width=width, height=height, halign=halign, valign=valign)
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
    color_num = [1,2,3,4],
    color_shades = fill(8, 4),
    show_icons::Bool = true,
    big_limits = (-2, 2, -4, 4),
    small_limits = nothing,  # set to a 4-tuple, or keep auto
    panel_rowgap = 8,
    panel_colgap = 10,
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
    cols = []
    for (i,palette) in enumerate(palettes)
        push!(cols, _make_colors(; palette=palette, n=nsets, color_num=color_num[i], color_shades=color_shades[i]))
    end

    # Layout inside the provided parent
    panel = GridLayout(parent)

    # Left: big axis spans two rows; Right: 2x2 grid of small axes
    ax_big = Axis(
        panel[1:2, 1],
        xlabel=L"\log_{10} \, \mu", xlabelsize=12,
        ylabel=L"\log_{10} \, \sigma^2", ylabelsize=12,
        xticklabelsize=12, yticklabelsize=12,
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

    # Create the 2x2 axes in the right grid
    axs_small = Axis[]
    grid_positions = [(2, 2), (2, 1), (1, 2), (1, 1)]
    for (i, (r, c)) in enumerate(grid_positions)
        i > nsets && break
        ax = Axis(
            right[r, c],
            xlabel=L"\log_{10} \, \mu", xlabelsize=12,
            ylabel=L"\log_{10} \, \sigma^2", ylabelsize=12,
            xticklabelsize=12, yticklabelsize=12,
        )
        if !isnothing(small_limits)
            xlims!(ax, small_limits[1], small_limits[2])
            ylims!(ax, small_limits[3], small_limits[4])
        end
        push!(axs_small, ax)
    end

    xtl = -10:0.1:10.0
    datasets = reverse(datasets)
    cols = reverse(cols)

    # Plot each dataset
    for (i, spec) in enumerate(datasets)
        i > length(axs_small) && break
        ax = axs_small[i]

        df = _load_tldf(spec.file)
        classes = unique(df.class)
        m_small_vec = []
        s_small_vec = []
        for (c,class) in enumerate(classes)
            sdf = df[df.class .== class, :]
            df_small = spec.occ_small > 0 ? sdf[sdf.occupancy .> spec.occ_small, :] : sdf
            df_small, m_small, s_small = _centered_logs(df_small)

            take = _ntget(spec, :take, :all)
            if take isa Integer
                n = min(take, length(m_small))
                idx = rand(1:length(m_small), n)
                m_small = m_small[idx]
                s_small = s_small[idx]
            end

            append!(m_small_vec, m_small)
            append!(s_small_vec, s_small)
    
            scatter!(
                ax, m_small, s_small;
                color=(:white, 1.0),
                strokecolor=cols[i][c],
                marker=markers[c],
                markersize=markersize_small,
                strokewidth=strokewidth,
            )
        end

        # same guide lines you used: y=2x and y=x, with the slight x-shift
        shift = _ntget(spec, :ref_shift, 0.0)
        lines!(ax, xtl .- minimum(m_small_vec) .+ shift, 2 .* (xtl .- minimum(m_small_vec)); linewidth=2, color=:black, linestyle=(:dash, :dense))
        lines!(ax, xtl .- minimum(m_small_vec) .+ shift, xtl; linewidth=2, color=:grey, linestyle=(:dash, :dense))

        # auto limits: your original scaling
        xlims!(ax, minimum(m_small_vec) * 1.5, maximum(m_small_vec) * 1.3)
        ylims!(ax, minimum(s_small_vec) * 1.5, maximum(s_small_vec) * 1.3)

        # icon overlay
        if show_icons && isfile(spec.icon)
            _add_icon!(right[grid_positions[i]...], spec.icon)
        end

        # big axis: occupancy filter + optional truncation
        for (c,class) in enumerate(classes)
            sdf = df[df.class .== class, :]
            df_big = spec.occ_big > 0 ? sdf[sdf.occupancy .> spec.occ_big, :] : sdf
            df_big, m_big, s_big = _centered_logs(df_big)
    
            take = _ntget(spec, :take, :all)
            if take isa Integer
                n = min(take, length(m_big))
                idx = rand(1:length(m_big), n)
                m_big = m_big[idx]
                s_big = s_big[idx]
            end
    
            scatter!(
                ax_big, m_big, s_big;
                color=(:white, 1.0),
                strokecolor=cols[i][c],
                marker=markers[c],
                markersize=markersize_big,
                strokewidth=strokewidth,
            )
        end
    end

    # Reference lines in all axes
    lines!(ax_big, xtl, 2 .* xtl; linewidth=2, color=:black, linestyle=(:dash, :dense))

    return (;
        panel,
        ax_big,
        axs_small,
        colors=cols,
    )
end

"""Backward-compatible wrapper that creates a standalone Figure."""
function plot_taylor(; TLDIR=Meris.DATADIR * "macro/taylor/", savefig=false, figname="tl.png", kwargs...)
    width = 1.9 * 246
    fig = Figure(; size=(width, width / 2), figure_padding=(2, 4, 2, 14))

    plot!(fig[1, 1]; TLDIR=TLDIR, kwargs...)

    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module
