#= Module to plot Figure 2B: the first Taylor's-law prediction bin panel across datasets =#
module Figure2B

using CairoMakie
using MakiePublication
using LaTeXStrings
using Colors
using DataFrames
using Statistics
using Meris

include(joinpath(@__DIR__, "plot-tl-prediction-bins.jl"))
using .TLPredictionBinPlotter

const MM_TO_PT = TLPredictionBinPlotter.MM_TO_PT
const NATURE_DOUBLE_WIDTH_PT = TLPredictionBinPlotter.NATURE_DOUBLE_WIDTH_PT
const NATURE_AXIS_LABEL_PT = TLPredictionBinPlotter.NATURE_AXIS_LABEL_PT
const NATURE_TICK_PT = TLPredictionBinPlotter.NATURE_TICK_PT
const FONT_SCALE = TLPredictionBinPlotter.FONT_SCALE
const LINE_WIDTH_SCALE = TLPredictionBinPlotter.LINE_WIDTH_SCALE
const MARKER_SIZE = TLPredictionBinPlotter.MARKER_SIZE

function _shades(base::Colorant, n::Int)
    hsl = convert(HSL, base)
    return [HSL(hsl.h, hsl.s, l) for l in range(0.25, 0.85, length=n)]
end

function _fig2_datasets(;
    TLDIR = Meris.DATADIR * "macro/taylor/",
    RESULTDIR = Meris.DATADIR * "macro/tl-prediction/",
)
    bases = [
        colorant"#1f77b4",
        colorant"#ff7f0e",
        colorant"#9467bd",
        colorant"#2ca02c",
        colorant"#d62728",
    ]

    palettes = [
        _shades(bases[1], 10),
        _shades(bases[2], 10),
        _shades(bases[3], 8),
        vcat(_shades(bases[4], 10)[1:7], _shades(bases[5], 8)),
    ]

    return [
        (;
            key="linguistic",
            title="Linguistic",
            prediction_file=joinpath(RESULTDIR, "linguistic.jld2"),
            palette=palettes[1],
            icon=joinpath(Meris.FIGDIR, "icons", "document.png"),
            icon_kw=(; width=Relative(0.24), height=Relative(0.24), halign=0.08, valign=0.92),
            take=1200,
        ),
        (;
            key="microbial",
            title="Microbial",
            prediction_file=joinpath(RESULTDIR, "microbial.jld2"),
            palette=palettes[2],
            icon=joinpath(Meris.FIGDIR, "icons", "bacteria.png"),
            icon_kw=(; width=Relative(0.24), height=Relative(0.24), halign=0.08, valign=0.92),
            take=1200,
        ),
        (;
            key="social",
            title="Social",
            prediction_file=joinpath(RESULTDIR, "social.jld2"),
            palette=palettes[3],
            icon=joinpath(Meris.FIGDIR, "icons", "socio-economic.png"),
            icon_kw=(; width=Relative(0.77 * 0.24), height=Relative(0.24), halign=0.08, valign=0.92),
            take=1600,
        ),
        (;
            key="biology",
            title="Biology",
            prediction_file=joinpath(RESULTDIR, "biology.jld2"),
            palette=palettes[4],
            icon=joinpath(Meris.FIGDIR, "icons", "eco.png"),
            icon_kw=(; width=Relative(0.24), height=Relative(0.24), halign=0.08, valign=0.92),
            take=1000,
        ),
    ]
end

function _prediction_curve(selected_rows::DataFrame, xmin, xmax)
    isempty(selected_rows.C_est) && return nothing

    logxs = range(xmin, xmax; length=300)
    xs = exp10.(logxs)
    C_est = median(selected_rows.C_est)
    return logxs, log10.(xs .+ C_est .* xs .^ 2)
end

function _plot_prediction_line!(ax, selected_rows::DataFrame, xmin, xmax)
    curve = _prediction_curve(selected_rows, xmin, xmax)
    isnothing(curve) && return ax
    logxs, logys = curve

    lines!(
        ax,
        logxs,
        logys;
        color=:black,
        linewidth=1.2 * LINE_WIDTH_SCALE,
        linestyle=(:dot, :dense),
    )

    return ax
end

function _plot_shaded_groups!(ax, groups, palette; markersize=MARKER_SIZE, strokewidth=0.65)
    for (j, group) in enumerate(groups)
        scatter!(
            ax,
            group.m,
            group.v;
            color=(:white, 1.0),
            strokecolor=palette[mod1(j, length(palette))],
            marker=group.marker,
            markersize=markersize,
            strokewidth=strokewidth,
        )
    end

    return ax
end

function _selected_bins_by_omega(component_bins::DataFrame, omega_log::Real; min_components::Int=10)
    bin_counts = combine(
        groupby(component_bins, [:class, :coeff_bin, :bin_center_log, :bin_center, :C_est, :C_fit]),
        nrow => :ncomponents,
    )
    filter!(:ncomponents => >=(min_components), bin_counts)
    filter!(:C_est => c -> isfinite(c) && c > 0, bin_counts)
    isempty(bin_counts.class) && return bin_counts

    selected = combine(groupby(bin_counts, :class)) do df
        distances = abs.(df.bin_center_log .- omega_log)
        candidates = df[distances .== minimum(distances), :]
        candidates[argmax(candidates.ncomponents), :]
    end
    sort!(selected, :class)

    return selected
end

function _data_xmax(groups; pad_fraction=0.04)
    isempty(groups) && return nothing
    xmax = maximum(maximum(group.m) for group in groups)
    xmin = minimum(minimum(group.m) for group in groups)
    return xmax + max((xmax - xmin) * pad_fraction, eps(Float64))
end

function _panel_limits(groups, selected_rows::DataFrame; xmin=-2.0, ypad_fraction=0.06)
    isempty(groups) && return nothing

    xmax = _data_xmax(groups)
    isnothing(xmax) && return nothing
    xmax = max(xmax, xmin + eps(Float64))

    ys = reduce(vcat, (group.v for group in groups))
    curve = _prediction_curve(selected_rows, xmin, xmax)
    if !isnothing(curve)
        _, curve_y = curve
        append!(ys, curve_y)
    end

    ymin, ymax = extrema(ys)
    ypad = max((ymax - ymin) * ypad_fraction, eps(Float64))
    return (xmin, xmax, ymin - ypad, ymax + ypad)
end

function _add_regime_background!(ax, panel_lims, threshold)
    xmin, xmax, ymin, ymax = panel_lims
    threshold >= xmax && return ax

    x0 = max(threshold, xmin)
    poly!(
        ax,
        Point2f[(x0, ymin), (xmax, ymin), (xmax, ymax), (x0, ymax)];
        color=(colorant"#d62728", 0.08),
        strokewidth=0,
    )

    return ax
end

function _add_linguistic_regime_labels!(ax, panel_lims, threshold; font_scale=1.0)
    xmin, xmax, _, _ = panel_lims
    fontsize = NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale * 0.75

    if threshold > xmin
        text!(
            ax,
            (xmin + threshold) / 2,
            2.3;
            text="linear",
            align=(:center, :center),
            fontsize=fontsize,
            color=:gray30,
        )
    end

    if threshold < xmax
        text!(
            ax,
            (threshold + xmax) / 2,
            -0.5;
            text="quadratic",
            align=(:center, :center),
            fontsize=fontsize,
            color=:gray30,
        )
    end

    return ax
end

function plot(;
    omega_log::Real=-0.5,
    xmin::Real=-2,
    datasets=_fig2_datasets(),
    ext::AbstractString="pdf",
    savefig::Bool=true,
    figname=Meris.FIGDIR * "fig2_B.$ext",
    font_scale::Float64=1.0,
    omega_font_scale::Float64=1.25,
    panel_colgap::Real=14,
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    fig = Figure(
        size=(1.18 * NATURE_DOUBLE_WIDTH_PT, 0.31 * NATURE_DOUBLE_WIDTH_PT),
        figure_padding=(6, 8, 6, 6),
    )

    for (i, spec) in enumerate(datasets)
        component_bins, _ = TLPredictionBinPlotter._load_bin_data(spec.prediction_file)
        selected_rows = _selected_bins_by_omega(component_bins, omega_log)
        groups = TLPredictionBinPlotter._bin_point_groups(component_bins, selected_rows; take=spec.take, seed_offset=100)

        ax = Axis(
            fig[1, i],
            xlabel=L"\log_{10}\,\mu",
            ylabel=L"\log_{10}\,\sigma^2",
            xlabelsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
            ylabelsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
            xticklabelsize=NATURE_TICK_PT * FONT_SCALE * font_scale,
            yticklabelsize=NATURE_TICK_PT * FONT_SCALE * font_scale,
            xgridvisible=false,
            ygridvisible=false,
        )

        panel_lims = _panel_limits(groups, selected_rows; xmin=Float64(xmin))
        threshold = -omega_log
        if !isnothing(panel_lims)
            _add_regime_background!(ax, panel_lims, threshold)
        end
        _plot_shaded_groups!(ax, groups, spec.palette; markersize=MARKER_SIZE)
        if !isnothing(panel_lims)
            _plot_prediction_line!(ax, selected_rows, panel_lims[1], panel_lims[2])
        end
        vlines!(
            ax,
            [threshold];
            color=:lightgrey,
            linewidth=0.8,
            linestyle=(:dash, :dense),
        )
        TLPredictionBinPlotter._add_icon!(fig[1, i], spec.icon; spec.icon_kw...)
        if !isnothing(panel_lims) && spec.key == "linguistic"
            _add_linguistic_regime_labels!(ax, panel_lims, threshold; font_scale=font_scale)
        end
        isnothing(panel_lims) || limits!(ax, panel_lims...)

        Label(
            fig[0, i],
            spec.title;
            fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
            tellwidth=false,
        )
    end

    Label(
        fig[-1, 1:4],
        latexstring("\\Omega \\sim 10^{", string(omega_log), "}");
        fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale * omega_font_scale,
        tellwidth=false,
    )

    rowgap!(fig.layout, 5)
    colgap!(fig.layout, panel_colgap)

    if savefig
        CairoMakie.save(figname, fig, pt_per_unit=1)
    end

    return fig
end

end # module Figure2B
