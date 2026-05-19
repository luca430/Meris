#= Composite Figure 2 plot, matching the two-tier sketch:
   A: Taylor's law decomposed by quadratic-coefficient bins.
   B: Dataset examples at one coefficient scale.
=#
module Fig2

using CairoMakie
using MakiePublication
using LaTeXStrings
using Colors
using DataFrames
using JLD2
using Random
using Statistics
using Meris

include(joinpath(@__DIR__, "plot-fig2-A.jl"))
using .Figure2A

include(joinpath(@__DIR__, "plot-fig2-B.jl"))
using .Figure2B

const TaylorPlotter = Figure2A.TaylorPlotter
const TLPredictionBinPlotter = Figure2B.TLPredictionBinPlotter

const NATURE_DOUBLE_WIDTH_PT = TLPredictionBinPlotter.NATURE_DOUBLE_WIDTH_PT
const NATURE_AXIS_LABEL_PT = TLPredictionBinPlotter.NATURE_AXIS_LABEL_PT
const NATURE_TICK_PT = TLPredictionBinPlotter.NATURE_TICK_PT
const NATURE_PANEL_LABEL_PT = TLPredictionBinPlotter.NATURE_PANEL_LABEL_PT
const FONT_SCALE = TLPredictionBinPlotter.FONT_SCALE
const LINE_WIDTH_SCALE = TLPredictionBinPlotter.LINE_WIDTH_SCALE
const MARKER_SIZE = TLPredictionBinPlotter.MARKER_SIZE

function _set_theme!()
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    set_theme!(MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true]))
end

function _fig2_datasets(;
    TLDIR=Meris.DATADIR * "macro/taylor/",
    RESULTDIR=Meris.DATADIR * "macro/tl-prediction/",
)
    datasets = Figure2B._fig2_datasets(; TLDIR=TLDIR, RESULTDIR=RESULTDIR)
    return [
        merge(spec, (; taylor_file=joinpath(TLDIR, "$(spec.key).jld2")))
        for spec in datasets
    ]
end

function _axis(parent; font_scale=1.0, xlabel=L"\log_{10}\,\mu", ylabel=L"\log_{10}\,\sigma^2")
    return Axis(
        parent;
        xlabel=xlabel,
        ylabel=ylabel,
        xlabelsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
        ylabelsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
        xticklabelsize=NATURE_TICK_PT * FONT_SCALE * font_scale,
        yticklabelsize=NATURE_TICK_PT * FONT_SCALE * font_scale,
        xgridvisible=false,
        ygridvisible=false,
    )
end

function _plot_taylor_law_panel!(
    parent,
    spec;
    font_scale=1.0,
    limits=(-1.8, 1.9, -1.8, 4.4),
    take=2200,
    markersize=MARKER_SIZE,
    strokewidth=0.65,
)
    ax = _axis(parent; font_scale=font_scale)
    df = TaylorPlotter._load_tldf(spec.taylor_file)
    if take isa Integer && nrow(df) > take
        Random.seed!(1234)
        df = df[randperm(nrow(df))[1:take], :]
    end

    markers = TLPredictionBinPlotter.BIN_MARKERS
    all_m = Float64[]
    all_v = Float64[]
    classes = unique(df.class)

    for (i, class) in enumerate(classes)
        sdf = df[df.class .== class, :]
        _, m, v = TaylorPlotter._log_mean_var(sdf; center=true)
        isempty(m) && continue

        append!(all_m, m)
        append!(all_v, v)
        scatter!(
            ax,
            m,
            v;
            color=(:white, 1.0),
            strokecolor=spec.palette[mod1(i, length(spec.palette))],
            marker=markers[mod1(i, length(markers))],
            markersize=markersize,
            strokewidth=strokewidth,
        )
    end

    if !isempty(all_m)
        xtl = range(limits[1], limits[2]; length=300)
        lines!(ax, xtl, 2 .* xtl; linewidth=1.3 * LINE_WIDTH_SCALE, color=:black, linestyle=(:dash, :dense))
        lines!(ax, xtl, xtl .- 0.35; linewidth=1.0 * LINE_WIDTH_SCALE, color=:gray55, linestyle=(:dash, :dense))
    end

    TLPredictionBinPlotter._add_icon!(parent, spec.icon; spec.icon_kw...)
    limits!(ax, limits...)
    return ax
end

function _plot_prediction_panel!(
    parent,
    spec,
    omega_log;
    font_scale=1.0,
    xmin=-1.6,
    markersize=MARKER_SIZE,
    show_labels=false,
)
    component_bins, _ = TLPredictionBinPlotter._load_bin_data(spec.prediction_file)
    selected_rows = Figure2B._selected_bins_by_omega(component_bins, omega_log)
    groups = TLPredictionBinPlotter._bin_point_groups(component_bins, selected_rows; take=spec.take, seed_offset=round(Int, 1000 + 100 * omega_log))
    ax = _axis(parent; font_scale=font_scale)
    panel_lims = Figure2B._panel_limits(groups, selected_rows; xmin=Float64(xmin))

    threshold = -omega_log
    if !isnothing(panel_lims)
        Figure2B._add_regime_background!(ax, panel_lims, threshold)
    end
    Figure2B._plot_shaded_groups!(ax, groups, spec.palette; markersize=markersize)
    if !isnothing(panel_lims)
        Figure2B._plot_prediction_line!(ax, selected_rows, panel_lims[1], panel_lims[2])
        vlines!(ax, [threshold]; color=:lightgrey, linewidth=0.8, linestyle=(:dash, :dense))
        show_labels && Figure2B._add_linguistic_regime_labels!(ax, panel_lims, threshold; font_scale=font_scale)
        limits!(ax, panel_lims...)
    end

    return ax
end

function _add_equation_symbol!(layout, row, col, text; font_scale=1.0)
    Label(
        layout[row, col],
        text;
        fontsize=NATURE_PANEL_LABEL_PT * FONT_SCALE * 1.55 * font_scale,
        tellwidth=true,
        tellheight=false,
    )
end

function _draw_timeline!(parent; font_scale=1.0)
    ax = Axis(parent; limits=(0, 1, 0, 1))
    hidedecorations!(ax)
    hidespines!(ax)

    lines!(ax, [0.06, 0.94], [0.50, 0.50]; color=:black, linewidth=1.2)
    scatter!(ax, [0.95], [0.50]; marker=:rtriangle, color=:black, markersize=15)
    for x in range(0.10, 0.88; length=12)
        lines!(ax, [x, x], [0.42, 0.58]; color=:black, linewidth=0.8)
    end

    for x in (0.24, 0.50, 0.76)
        lines!(ax, [x - 0.13, x, x + 0.13], [0.96, 0.50, 0.96]; color=:gray70, linewidth=0.6)
    end

    text!(ax, 0.975, 0.54; text=L"\Omega", fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * 1.25 * font_scale, align=(:left, :center), color=:black)
    text!(ax, 0.945, 0.18; text="quadratic\ncoefficient", fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * 0.72 * font_scale, align=(:left, :center), color=:black)

    return ax
end

function _plot_lower_panels!(parent, datasets, omega_log; font_scale=1.0, panel_colgap=14)
    lower = GridLayout(parent)
    Box(lower[1:3, 1:4]; color=(:white, 0), strokecolor=:gray70, strokewidth=0.7)

    Label(
        lower[1, 1:4],
        latexstring("\\Omega \\sim 10^{", string(omega_log), "}");
        fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale * 1.25,
        tellwidth=false,
    )
    Label(
        lower[1, 1, TopLeft()],
        "B";
        fontsize=NATURE_PANEL_LABEL_PT * FONT_SCALE * font_scale,
        font=:bold,
        padding=(4, 0, 0, 0),
        halign=:left,
        valign=:top,
        tellwidth=false,
        tellheight=false,
    )

    for (i, spec) in enumerate(datasets)
        Label(
            lower[2, i],
            spec.title;
            fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
            tellwidth=false,
        )
        _plot_prediction_panel!(
            lower[3, i],
            spec,
            omega_log;
            font_scale=font_scale,
            xmin=-2.0,
            markersize=MARKER_SIZE,
            show_labels=spec.key == "linguistic",
        )
        TLPredictionBinPlotter._add_icon!(lower[3, i], spec.icon; spec.icon_kw...)
    end

    rowgap!(lower, 1, 1)
    rowgap!(lower, 2, 3)
    colgap!(lower, panel_colgap)
    return lower
end

function plot(;
    ext::AbstractString="pdf",
    savefig::Bool=true,
    figname=Meris.FIGDIR * "fig2.$ext",
    datasets=_fig2_datasets(),
    upper_omegas=(-1, 0, 1),
    lower_omega=-0.5,
    font_scale::Float64=1.0,
)
    _set_theme!()

    fig = Figure(
        size=(1.55 * NATURE_DOUBLE_WIDTH_PT, 0.83 * NATURE_DOUBLE_WIDTH_PT),
        figure_padding=(8, 10, 6, 8),
    )

    linguistic = datasets[findfirst(spec -> spec.key == "linguistic", datasets)]
    top = GridLayout(fig[1, 1])

    Label(
        top[1, 1, TopLeft()],
        "A";
        fontsize=NATURE_PANEL_LABEL_PT * FONT_SCALE * font_scale,
        font=:bold,
        padding=(0, 0, 0, 0),
        halign=:left,
        valign=:top,
        tellwidth=false,
        tellheight=false,
    )
    Label(
        top[1, 1],
        "Taylor's law";
        fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
        tellwidth=false,
    )
    _plot_taylor_law_panel!(top[2, 1], linguistic; font_scale=font_scale)

    _add_equation_symbol!(top, 2, 2, "="; font_scale=font_scale)
    symbols = ["+", "+"]
    for (j, omega_log) in enumerate(upper_omegas)
        axis_col = 2 * j + 1
        Label(
            top[1, axis_col],
            latexstring("\\Omega \\sim 10^{", string(omega_log), "}");
            fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
            tellwidth=false,
        )
        Box(top[2, axis_col]; color=(:white, 0), strokecolor=:gray70, strokewidth=0.7)
        _plot_prediction_panel!(
            top[2, axis_col],
            linguistic,
            omega_log;
            font_scale=font_scale,
            xmin=-1.6,
            markersize=MARKER_SIZE,
            show_labels=omega_log == 0,
        )
        if j <= length(symbols)
            _add_equation_symbol!(top, 2, axis_col + 1, symbols[j]; font_scale=font_scale)
        end
    end

    _add_equation_symbol!(top, 2, 8, "+"; font_scale=font_scale)
    Label(
        top[2, 9],
        "...";
        fontsize=NATURE_PANEL_LABEL_PT * FONT_SCALE * 1.35 * font_scale,
        tellwidth=true,
        tellheight=false,
    )

    for c in (2, 4, 6, 8)
        colsize!(top, c, Fixed(22))
    end
    colsize!(top, 9, Fixed(32))
    rowgap!(top, 1, 2)
    colgap!(top, 7)

    _draw_timeline!(fig[2, 1]; font_scale=font_scale)
    _plot_lower_panels!(fig[3, 1], datasets, lower_omega; font_scale=font_scale)

    rowsize!(fig.layout, 1, Relative(0.44))
    rowsize!(fig.layout, 2, Relative(0.13))
    rowsize!(fig.layout, 3, Relative(0.43))
    rowgap!(fig.layout, 4)

    if savefig
        CairoMakie.save(figname, fig, pt_per_unit=1)
    end

    return fig
end

end # module Fig2
