#= Module to plot held-out mean/variance panels for Taylor's-law prediction bins =#
#/ Start module
module TLPredictionBinPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings
using JLD2
using DataFrames
using Colors
using FileIO, ImageTransformations
using Random
using Statistics

#/ Modules
using Meris

const ICONDIR = Meris.FIGDIR .* "icons"
const MM_TO_PT = 72.0 / 25.4
const NATURE_SINGLE_WIDTH_PT = 89.0 * MM_TO_PT
const NATURE_DOUBLE_WIDTH_PT = 183.0 * MM_TO_PT
const NATURE_AXIS_LABEL_PT = 7
const NATURE_TICK_PT = 6
const NATURE_TEXT_PT = 6
const NATURE_PANEL_LABEL_PT = 8
const FONT_SCALE = 1.75
const LINE_WIDTH_SCALE = 2.0
const MARKER_SIZE = 6.2

function _shades(base::Colorant, n::Int)
    hsl = convert(HSL, base)
    return [HSL(hsl.h, hsl.s, l) for l in range(0.25, 0.85, length=n)]
end

function _load_bin_data(path::AbstractString)
    d = JLD2.load(path)
    return d["component_bins"], d["selected_bins"]
end

function _load_tldf(path::AbstractString)
    return JLD2.load(path)["tldf"]
end

function _centered_logs(df; mean_col=:omeanfrequency, var_col=:ovarfrequency)
    d = df[df[!, var_col] .> 0.0, :]
    m = log10.(d[!, mean_col])
    v = log10.(d[!, var_col])
    m .-= mean(m)
    v .-= mean(v)
    return d, m, v
end

function _default_datasets(;
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

    return [
        (;
            key="linguistic",
            title="Linguistic",
            taylor_file=joinpath(TLDIR, "linguistic.jld2"),
            prediction_file=joinpath(RESULTDIR, "linguistic.jld2"),
            palette=_shades(bases[1], 10),
            icon=joinpath(ICONDIR, "document.png"),
            icon_kw=(; width=Relative(0.24), height=Relative(0.24), halign=0.08, valign=0.92),
            take=1200,
        ),
        (;
            key="microbial",
            title="Microbial",
            taylor_file=joinpath(TLDIR, "microbial.jld2"),
            prediction_file=joinpath(RESULTDIR, "microbial.jld2"),
            palette=_shades(bases[2], 10),
            icon=joinpath(ICONDIR, "bacteria.png"),
            icon_kw=(; width=Relative(0.24), height=Relative(0.24), halign=0.08, valign=0.92),
            take=1200,
        ),
        (;
            key="social",
            title="Social",
            taylor_file=joinpath(TLDIR, "social.jld2"),
            prediction_file=joinpath(RESULTDIR, "social.jld2"),
            palette=_shades(bases[3], 8),
            icon=joinpath(ICONDIR, "socio-economic.png"),
            icon_kw=(; width=Relative(0.77 * 0.24), height=Relative(0.24), halign=0.08, valign=0.92),
            take=1600,
        ),
        (;
            key="biology",
            title="Biology",
            taylor_file=joinpath(TLDIR, "biology.jld2"),
            prediction_file=joinpath(RESULTDIR, "biology.jld2"),
            palette=vcat(_shades(bases[5], 8)[1:4], _shades(bases[4], 10)[1:7]),
            icon=joinpath(ICONDIR, "eco.png"),
            icon_kw=(; width=Relative(0.24), height=Relative(0.24), halign=0.08, valign=0.92),
            take=1000,
        ),
    ]
end

function _add_icon!(parent_cell, icon_path; width=Relative(0.24), height=Relative(0.24), halign=0.08, valign=0.92)
    isfile(icon_path) || return nothing

    axicon = Axis(
        parent_cell;
        width=width,
        height=height,
        halign=halign,
        valign=valign,
        tellwidth=false,
        tellheight=false,
    )

    icon = FileIO.load(icon_path)
    icon_small = imresize(icon, (320, 320))
    image!(axicon, rotr90(icon_small))
    hidedecorations!(axicon)
    hidespines!(axicon)

    return axicon
end

function _plot_taylor_small!(ax, df::DataFrame, palette; take=:all, markersize=MARKER_SIZE, strokewidth=0.65)
    if take isa Integer && nrow(df) > take
        Random.seed!(1234)
        df = df[randperm(nrow(df))[1:take], :]
    end

    markers = [:circle, :rect, :diamond, :cross, :x, :utriangle, :dtriangle, :star4, :star6, :pentagon, :hexagon, :octagon]
    classes = unique(df.class)
    all_m = Float64[]
    all_v = Float64[]

    for (i, class) in enumerate(classes)
        sdf = df[df.class .== class, :]
        _, m, v = _centered_logs(sdf)
        isempty(m) && continue

        append!(all_m, m)
        append!(all_v, v)
        scatter!(
            ax,
            m,
            v;
            color=(:white, 1.0),
            strokecolor=palette[mod1(i, length(palette))],
            marker=markers[mod1(i, length(markers))],
            markersize=markersize,
            strokewidth=strokewidth,
        )
    end

    if !isempty(all_m)
        xtl = range(minimum(all_m), maximum(all_m); length=200)
        x0 = minimum(all_m)
        y0 = minimum(all_v)
        lines!(ax, xtl, 2 .* (xtl .- x0) .+ y0; linewidth=1.2 * LINE_WIDTH_SCALE, color=:black, linestyle=(:dash, :dense))
        lines!(ax, xtl, (xtl .- x0) .+ y0; linewidth=1.0 * LINE_WIDTH_SCALE, color=:grey, linestyle=(:dash, :dense))
    end

    return ax
end

function _plot_taylor_component_bins!(
    ax,
    component_bins::DataFrame,
    palette;
    take=:all,
    markersize=MARKER_SIZE,
    strokewidth=0.65,
    black_yshift=0.0,
    gray_yshift=0.0,
)
    markers = [:circle, :rect, :diamond, :cross, :x, :utriangle, :dtriangle, :star4, :star6, :pentagon, :hexagon, :octagon]
    classes = unique(component_bins.class)
    all_m = Float64[]
    all_v = Float64[]
    max_points_per_class = take isa Integer ?
        max(1, ceil(Int, take / length(classes))) :
        typemax(Int)

    for (i, class) in enumerate(classes)
        d = component_bins[component_bins.class .== class, :]
        filter!([:mean, :var] => (m, v) -> isfinite(m) && isfinite(v) && m > 0 && v > 0, d)
        isempty(d.mean) && continue

        if nrow(d) > max_points_per_class
            Random.seed!(1234 + i)
            d = d[randperm(nrow(d))[1:max_points_per_class], :]
        end

        m = log10.(d.mean)
        v = log10.(d.var)

        append!(all_m, m)
        append!(all_v, v)
        scatter!(
            ax,
            m,
            v;
            color=(:white, 1.0),
            strokecolor=palette[mod1(i, length(palette))],
            marker=markers[mod1(i, length(markers))],
            markersize=markersize,
            strokewidth=strokewidth,
        )
    end

    if !isempty(all_m)
        xtl = range(minimum(all_m), maximum(all_m); length=200)
        x0 = minimum(all_m)
        y0 = minimum(all_v)
        lines!(ax, xtl, 2 .* (xtl .- x0) .+ y0 .+ black_yshift; linewidth=1.2 * LINE_WIDTH_SCALE, color=:black, linestyle=(:dash, :dense))
        lines!(ax, xtl, (xtl .- x0) .+ y0 .+ gray_yshift; linewidth=1.0 * LINE_WIDTH_SCALE, color=:grey, linestyle=(:dash, :dense))
    end

    return ax
end

function _selected_bin(selected_bins::DataFrame, coeff_order::Int)
    d = selected_bins[selected_bins.coeff_order .== coeff_order, :]
    isempty(d.coeff_order) && return nothing
    return d[argmax(d.ncomponents), :]
end

function _selected_bins_by_class(component_bins::DataFrame, coeff_order::Int; min_components::Int=10)
    bin_counts = combine(
        groupby(component_bins, [:class, :coeff_order, :coeff_bin, :bin_center_log, :bin_center, :C_est, :C_fit]),
        nrow => :ncomponents,
    )
    filter!(:coeff_order => ==(coeff_order), bin_counts)
    filter!(:ncomponents => >=(min_components), bin_counts)
    filter!(:C_est => c -> isfinite(c) && c > 0, bin_counts)
    isempty(bin_counts.class) && return bin_counts

    selected = combine(groupby(bin_counts, :class)) do df
        df[argmax(df.ncomponents), :]
    end
    sort!(selected, :class)

    return selected
end

function _plot_bin_panel!(ax, component_bins::DataFrame, selected_row; markersize=MARKER_SIZE, strokewidth=0.65)
    isnothing(selected_row) && return ax

    d = component_bins[
        (component_bins.class .== selected_row.class) .&
        (component_bins.coeff_bin .== selected_row.coeff_bin),
        :,
    ]
    filter!([:mean, :var] => (m, v) -> isfinite(m) && isfinite(v) && m > 0 && v > 0, d)
    isempty(d.mean) && return ax

    scatter!(
        ax,
        log10.(d.mean),
        log10.(d.var);
        color=(:white, 1.0),
        strokecolor=:black,
        markersize=markersize,
        strokewidth=strokewidth,
    )

    logxmin, logxmax = extrema(log10.(d.mean))
    logxs = range(logxmin, logxmax; length=300)
    xs = exp10.(logxs)

    lines!(
        ax,
        logxs,
        log10.(xs .+ selected_row.C_est .* xs .^ 2);
        color=:black,
        linewidth=1.2 * LINE_WIDTH_SCALE,
        linestyle=(:dot, :dense),
    )

    return ax
end

function _plot_bin_panel_all_classes!(
    ax,
    component_bins::DataFrame,
    selected_rows::DataFrame,
    palette;
    take=:all,
    markersize=MARKER_SIZE,
    strokewidth=0.65,
)
    isempty(selected_rows.class) && return ax

    all_means = Float64[]
    all_vars = Float64[]
    markers = [:circle, :rect, :diamond, :cross, :x, :utriangle, :dtriangle, :star4, :star6, :pentagon, :hexagon, :octagon]
    max_points_per_class = take isa Integer ?
        max(1, ceil(Int, take / nrow(selected_rows))) :
        typemax(Int)

    for (i, row) in enumerate(eachrow(selected_rows))
        d = component_bins[
            (component_bins.class .== row.class) .&
            (component_bins.coeff_bin .== row.coeff_bin),
            :,
        ]
        filter!([:mean, :var] => (m, v) -> isfinite(m) && isfinite(v) && m > 0 && v > 0, d)
        isempty(d.mean) && continue

        if nrow(d) > max_points_per_class
            Random.seed!(1234 + i)
            d = d[randperm(nrow(d))[1:max_points_per_class], :]
        end

        append!(all_means, d.mean)
        append!(all_vars, d.var)
    scatter!(
        ax,
        log10.(d.mean),
        log10.(d.var);
            color=(:white, 1.0),
            strokecolor=palette[mod1(i, length(palette))],
            marker=markers[mod1(i, length(markers))],
            markersize=markersize,
            strokewidth=strokewidth,
        )
    end

    isempty(all_means) && return ax

    logxmin, logxmax = extrema(log10.(all_means))
    logxs = range(logxmin, logxmax; length=300)
    xs = exp10.(logxs)
    C_est = median(selected_rows.C_est)

    lines!(
        ax,
        logxs,
        log10.(xs .+ C_est .* xs .^ 2);
        color=:black,
        linewidth=1.2 * LINE_WIDTH_SCALE,
        linestyle=(:dot, :dense),
    )

    return ax
end

function plot_test_bins(
    component_bins::DataFrame,
    selected_bins::DataFrame;
    savefig::Bool=true,
    filename=Meris.FIGDIR * "gutenberg-gtex-test-bin-mean-var.pdf",
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    n = nrow(selected_bins)
    ncols = min(3, n)
    nrows = ceil(Int, n / ncols)

    fig = Figure(
        size = (ncols * 0.82 * NATURE_SINGLE_WIDTH_PT, nrows * 0.72 * NATURE_SINGLE_WIDTH_PT),
        figure_padding = (4, 8, 4, 4),
    )

    for (i, row) in enumerate(eachrow(selected_bins))
        r = div(i - 1, ncols) + 1
        c = mod(i - 1, ncols) + 1

        ax = Axis(
            fig[r, c],
            xlabel = L"\log_{10}\,\mu",
            ylabel = L"\log_{10}\,\sigma^2",
            xlabelsize = NATURE_AXIS_LABEL_PT * FONT_SCALE,
            ylabelsize = NATURE_AXIS_LABEL_PT * FONT_SCALE,
            xticklabelsize = NATURE_TICK_PT * FONT_SCALE,
            yticklabelsize = NATURE_TICK_PT * FONT_SCALE,
            xgridvisible = false,
            ygridvisible = false,
        )

        d = component_bins[
            (component_bins.class .== row.class) .&
            (component_bins.coeff_bin .== row.coeff_bin),
            :,
        ]
        filter!([:mean, :var] => (m, v) -> isfinite(m) && isfinite(v) && m > 0 && v > 0, d)
        isempty(d.mean) && continue

        scatter!(
            ax,
            log10.(d.mean),
            log10.(d.var);
            color = (:white, 1.0),
            strokecolor = :black,
            markersize = MARKER_SIZE,
            strokewidth = 0.7,
        )

        logxmin, logxmax = extrema(log10.(d.mean))
        logxs = range(logxmin, logxmax; length=300)
        xs = exp10.(logxs)

        lines!(
            ax,
            logxs,
            log10.(xs .+ row.C_est .* xs .^ 2);
            color = :black,
            linewidth = 1.4 * LINE_WIDTH_SCALE,
            linestyle = (:dot, :dense),
        )
    end

    if savefig
        CairoMakie.save(filename, fig, pt_per_unit=1)
    end

    return fig
end

function plot_grid(;
    datasets=_default_datasets(),
    coeff_orders=(-1, 0, 1),
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    font_scale::Float64=1.0,
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    nrows = length(datasets)
    axis_cols = [1, 3, 5, 7]
    symbol_cols = [2, 4, 6]
    symbols = ["=", "+", "+"]

    fig = Figure(
        size=(1.30 * NATURE_DOUBLE_WIDTH_PT, 1.02 * NATURE_DOUBLE_WIDTH_PT),
        figure_padding=(6, 10, 6, 6),
    )

    for (r, spec) in enumerate(datasets)
        component_bins, selected_bins = _load_bin_data(spec.prediction_file)

        ax_tl = Axis(
            fig[r, axis_cols[1]],
            xlabel = r == nrows ? L"\log_{10}\,\mu" : "",
            ylabel = L"\log_{10}\,\sigma^2",
            xlabelsize = NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
            ylabelsize = NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
            xticklabelsize = NATURE_TICK_PT * FONT_SCALE * font_scale,
            yticklabelsize = NATURE_TICK_PT * FONT_SCALE * font_scale,
            xgridvisible=false,
            ygridvisible=false,
        )
        _plot_taylor_component_bins!(ax_tl, component_bins, spec.palette; take=spec.take, black_yshift=0.4, gray_yshift=-0.55)
        _add_icon!(fig[r, axis_cols[1]], spec.icon; spec.icon_kw...)

        for (j, coeff_order) in enumerate(coeff_orders)
            c = axis_cols[j + 1]
            ax = Axis(
                fig[r, c],
                xlabel = r == nrows ? L"\log_{10}\,\mu" : "",
                ylabel = j == 1 ? L"\log_{10}\,\sigma^2" : "",
                xlabelsize = NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
                ylabelsize = NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
                xticklabelsize = NATURE_TICK_PT * FONT_SCALE * font_scale,
                yticklabelsize = NATURE_TICK_PT * FONT_SCALE * font_scale,
                xgridvisible=false,
                ygridvisible=false,
            )
            selected_rows = _selected_bins_by_class(component_bins, coeff_order)
            _plot_bin_panel_all_classes!(ax, component_bins, selected_rows, spec.palette; take=spec.take)

            if r == 1
                Label(
                    fig[0, c],
                    latexstring("\\Omega \\sim 10^{", string(coeff_order), "}");
                    fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
                    tellwidth=false,
                )
            end
        end

        for (c, sym) in zip(symbol_cols, symbols)
            Label(
                fig[r, c],
                sym;
                fontsize=NATURE_PANEL_LABEL_PT * 2.1 * font_scale,
                tellwidth=true,
                tellheight=false,
            )
        end
    end

    Label(
        fig[0, axis_cols[1]],
        "Taylor's law";
        fontsize=NATURE_AXIS_LABEL_PT * FONT_SCALE * font_scale,
        tellwidth=false,
    )

    for c in symbol_cols
        colsize!(fig.layout, c, Fixed(18))
    end

    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 8)

    if savefig
        outfile = isnothing(figname) ? (Meris.FIGDIR * "tl-prediction-bin-grid.$ext") : figname
        CairoMakie.save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

function plot_category(;
    RESULTDIR = Meris.DATADIR * "macro/tl-prediction/",
    category = "linguistic",
    FILENAME = nothing,
    ext = "pdf",
    savefig::Bool = true,
    figname = nothing,
    kwargs...
)
    file = isnothing(FILENAME) ? "$(category).jld2" : FILENAME
    component_bins, selected_bins = _load_bin_data(RESULTDIR * file)
    outfile = isnothing(figname) ? (Meris.FIGDIR * "tl-prediction-bins-$(category).$ext") : figname
    return plot_test_bins(component_bins, selected_bins; savefig=savefig, filename=outfile, kwargs...)
end

function plot(; kwargs...)
    return plot_grid(; kwargs...)
end

end # module TLPredictionBinPlotter
#/ End module
