#= Module to plot Taylor's-law prediction summaries across Figure 2 datasets =#
#/ Start module
module TLPredictionPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings
using JLD2
using DataFrames
using Colors
using Statistics
using FileIO

#/ Modules
using Meris

const ICONDIR = Meris.FIGDIR .* "icons"
const MM_TO_PT = 72.0 / 25.4
const NATURE_SINGLE_WIDTH_PT = 89.0 * MM_TO_PT
const NATURE_DOUBLE_WIDTH_PT = 183.0 * MM_TO_PT
const NATURE_AXIS_LABEL_PT = 7
const NATURE_TICK_PT = 6
const NATURE_TEXT_PT = 6

const TL_MARKERS = [:circle, :rect, :diamond, :cross, :x, :utriangle, :dtriangle, :star4, :star6, :pentagon, :hexagon, :octagon]

#################
### INTERNALS ###

function _shades(base::Colorant, n::Int)
    hsl = convert(HSL, base)
    return [HSL(hsl.h, hsl.s, l) for l in range(0.25, 0.85, length=n)]
end

function _load_prediction_data(path::AbstractString)
    d = JLD2.load(path)

    if haskey(d, "plot_summaries")
        return d["plot_summaries"], d["component_bins"], :fit_ab
    elseif haskey(d, "fit_summaries")
        return d["fit_summaries"], d["component_bins"], :fit_ab
    end

    return d["test_summaries"], d["component_bins"], :legacy
end

function _valid_bins(component_bins::DataFrame, class; min_components::Int)
    d = component_bins[component_bins.class .== String(class), :]

    isempty(d.class) && return Set{Int}()

    bin_counts = combine(groupby(d, :coeff_bin), nrow => :ncomponents)
    filter!(:ncomponents => >=(min_components), bin_counts)

    return Set(Int.(bin_counts.coeff_bin))
end

function _default_datasets(; RESULTDIR = Meris.DATADIR * "macro/tl-prediction/")
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
        vcat(_shades(bases[5], 8)[1:4], _shades(bases[4], 10)[1:7]),
    ]

    return [
        (;
            key="linguistic",
            title="Linguistic",
            prediction_file=joinpath(RESULTDIR, "linguistic.jld2"),
            palette=palettes[1],
            icon=joinpath(ICONDIR, "document.png"),
            icon_kw=(; width=Relative(0.16), height=Relative(0.16), halign=0.18, valign=0.92),
        ),
        (;
            key="microbial",
            title="Microbial",
            prediction_file=joinpath(RESULTDIR, "microbial.jld2"),
            palette=palettes[2],
            icon=joinpath(ICONDIR, "bacteria.png"),
            icon_kw=(; width=Relative(0.16), height=Relative(0.16), halign=0.18, valign=0.92),
        ),
        (;
            key="social",
            title="Social",
            prediction_file=joinpath(RESULTDIR, "social.jld2"),
            palette=palettes[3],
            icon=joinpath(ICONDIR, "socio-economic.png"),
            icon_kw=(; width=Relative(0.77 * 0.16), height=Relative(0.16), halign=0.18, valign=0.92),
        ),
        (;
            key="biology",
            title="Biology",
            prediction_file=joinpath(RESULTDIR, "biology.jld2"),
            palette=palettes[4],
            icon=joinpath(ICONDIR, "eco.png"),
            icon_kw=(; width=Relative(0.20), height=Relative(0.20), halign=0.18, valign=0.92),
        ),
    ]
end

function _select_datasets(dataset, all_datasets)
    if dataset isa AbstractString || dataset isa Symbol
        key = lowercase(String(dataset))
        key == "all" && return all_datasets
        wanted = [key]
    else
        wanted = lowercase.(String.(dataset))
    end

    available = [lowercase(String(spec.key)) for spec in all_datasets]
    selected = [spec for spec in all_datasets if lowercase(String(spec.key)) in wanted]
    missing = setdiff(wanted, available)
    isempty(missing) || error("Unknown dataset(s): $(join(missing, ", ")). Available datasets: all, $(join(available, ", "))")

    return selected
end

function _class_name(class)
    s = String(class)
    s = replace(s, r"^gutenberg-"i => "")
    s = replace(s, r"^guten-"i => "")
    s = replace(s, r"^gen-"i => "")
    return s
end

function _ntget(nt::NamedTuple, field::Symbol, default)
    return hasproperty(nt, field) ? getproperty(nt, field) : default
end

function _add_icon!(parent_cell, icon_path;
    width=Relative(0.20), height=Relative(0.20),
    halign=0.18, valign=0.95
)
    axicon = Axis(
        parent_cell;
        width=width, height=height,
        halign=halign, valign=valign,
        tellwidth=false, tellheight=false,
        aspect=DataAspect(),
    )

    icon = FileIO.load(icon_path)
    image!(axicon, rotr90(icon))
    hidedecorations!(axicon)
    hidespines!(axicon)

    return axicon
end

function _plot_records(datasets; min_components::Int=11)
    records = NamedTuple[]

    for spec in datasets
        summaries, component_bins, summary_kind = _load_prediction_data(spec.prediction_file)
        classes = if String(spec.key) == "biology"
            sort(
                collect(keys(summaries));
                by=class -> (startswith(String(class), "gen-") ? 0 : 1, string(class)),
            )
        else
            sort(collect(keys(summaries)); by=string)
        end

        for (i, class) in enumerate(classes)
            valid_bins = _valid_bins(component_bins, class; min_components=min_components)
            df = copy(summaries[class])
            df = df[in.(df.coeff_bin, Ref(valid_bins)), :]

            if summary_kind == :fit_ab
                df.C_est = df.C_fit_A
                df.C_fit = df.C_fit_B
                df.C_est_err = df.C_fit_A_err
                df.C_fit_err = df.C_fit_B_err
            end

            push!(
                records,
                (;
                    dataset=String(spec.title),
                    dataset_key=String(spec.key),
                    icon=_ntget(spec, :icon, nothing),
                    icon_kw=_ntget(spec, :icon_kw, (;)),
                    class=class,
                    df=df,
                    color=spec.palette[mod1(i, length(spec.palette))],
                    marker=TL_MARKERS[mod1(i, length(TL_MARKERS))],
                )
            )
        end
    end

    return records
end

function _clean_record_data(record)
    d = copy(record.df)
    filter!(
        [:C_est, :C_fit] =>
            (x, y) -> all(isfinite, (x, y)) && x > 0 && y > 0,
        d
    )
    return d
end

function _fit_slope(records)
    slopes = Float64[]

    for record in records
        d = _clean_record_data(record)
        nrow(d) == 0 && continue

        x = Float64.(d.C_est)
        y = Float64.(d.C_fit)
        denom = sum(abs2, x)
        denom > 0 || continue

        push!(slopes, sum(x .* y) / denom)
    end

    isempty(slopes) && return nothing

    return (;
        a=mean(slopes),
        err=length(slopes) > 1 ? std(slopes) / sqrt(length(slopes)) : 0.0,
        n=length(slopes),
    )
end

function _record_groups_by_dataset(records)
    order = String[]
    groups = Dict{String, Vector{NamedTuple}}()

    for record in records
        key = record.dataset_key
        if !haskey(groups, key)
            groups[key] = NamedTuple[]
            push!(order, key)
        end
        push!(groups[key], record)
    end

    return [(key, groups[key]) for key in order]
end

function _global_limits(records; logscale::Bool=true)
    values = Float64[]

    for record in records
        d = _clean_record_data(record)
        nrow(d) == 0 && continue
        append!(values, d.C_est)
        append!(values, d.C_fit)
    end

    isempty(values) && return nothing

    lo = minimum(values)
    hi = maximum(values)
    if lo == hi
        lo, hi = logscale ? (lo / 10, hi * 10) : (lo - 1, hi + 1)
    end

    return (lo, hi)
end

function _plot_dataset_records!(
    ax,
    records;
    errorbars::Bool=false,
    xerr_col::Symbol=:C_est_err,
    yerr_col::Symbol=:C_fit_err,
    logscale::Bool=true,
    markersize::Int=6,
    strokewidth::Float64=0.7,
    font_scale::Float64=2.0,
    limits=nothing,
)
    any_points = false

    for record in records
        d = _clean_record_data(record)

        if nrow(d) == 0
            @warn "Skipping class $(record.class): no valid points"
            continue
        end

        any_points = true
        x = d.C_est
        y = d.C_fit

        if errorbars
            d_err = copy(d)
            filter!(
                [:C_est, :C_fit, xerr_col, yerr_col] =>
                    (xv, yv, xe, ye) ->
                        all(isfinite, (xv, yv, xe, ye)) &&
                        xe >= 0 &&
                        ye >= 0 &&
                        xv - xe > 0 &&
                        yv - ye > 0,
                d_err
            )

            if nrow(d_err) > 0
                errorbars!(
                    ax,
                    d_err.C_est,
                    d_err.C_fit,
                    d_err[:, yerr_col];
                    color = record.color,
                    whiskerwidth = 4,
                    linewidth = 0.7,
                )

                errorbars!(
                    ax,
                    d_err.C_est,
                    d_err.C_fit,
                    d_err[:, xerr_col];
                    direction = :x,
                    color = record.color,
                    whiskerwidth = 4,
                    linewidth = 0.7,
                )
            end
        end

        scatter!(
            ax,
            x,
            y;
            marker = record.marker,
            color = (:white, 1.0),
            strokecolor = record.color,
            markersize = markersize,
            strokewidth = strokewidth,
        )
    end

    fit = _fit_slope(records)

    if any_points && !isnothing(fit)
        lo, hi = isnothing(limits) ? _global_limits(records; logscale=logscale) : limits
        xs = logscale ? exp10.(range(log10(lo), log10(hi); length=300)) :
                        range(lo, hi; length=300)

        lines!(
            ax,
            xs,
            fit.a .* xs;
            linestyle = (:dash, :dense),
            color = :black,
            linewidth = 2.2,
        )

        text!(
            ax,
            0.95,
            0.10;
            text = L"a = %$(round(fit.a; digits=2)) \pm %$(round(fit.err; digits=2))",
            space = :relative,
            align = (:right, :bottom),
            fontsize = NATURE_TEXT_PT * font_scale,
            color = :black,
        )
    end

    return (; ax, fit, any_points)
end

#################
### FUNCTIONS ###

function plot_C_est_vs_C_fit(
    records;
    errorbars::Bool = false,
    xerr_col::Symbol = :C_est_err,
    yerr_col::Symbol = :C_fit_err,
    logscale::Bool = true,
    markersize::Int = 8,
    strokewidth::Float64 = 0.7,
    font_scale::Float64 = 2.4,
    show_icons::Bool = true,
    xlimits = (nothing, 1e2),
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    dataset_groups = _record_groups_by_dataset(records)
    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]
    nplots = length(dataset_groups)
    nrows = cld(max(nplots, 1), 2)
    limits = _global_limits(records; logscale=logscale)
    plot_xlimits = if isnothing(limits)
        nothing
    else
        lo, hi = limits
        (
            isnothing(xlimits[1]) ? lo : xlimits[1],
            isnothing(xlimits[2]) ? hi : xlimits[2],
        )
    end

    fig = Figure(
        size = (NATURE_DOUBLE_WIDTH_PT, 0.78 * NATURE_DOUBLE_WIDTH_PT),
        figure_padding = (4, 6, 4, 4)
    )

    axes = Axis[]

    for (i, (_, group_records)) in enumerate(dataset_groups)
        i <= length(positions) || break
        row, col = positions[i]

        ax = Axis(
            fig[row, col],
            xlabel = row == nrows ? L"\Omega_{\mathrm{fit},A}" : "",
            ylabel = col == 1 ? L"\Omega_{\mathrm{fit},B}" : "",
            xscale = logscale ? log10 : identity,
            yscale = logscale ? log10 : identity,
            xlabelsize = NATURE_AXIS_LABEL_PT * font_scale,
            ylabelsize = NATURE_AXIS_LABEL_PT * font_scale,
            xticklabelsize = NATURE_TICK_PT * font_scale,
            yticklabelsize = NATURE_TICK_PT * font_scale,
            aspect = AxisAspect(1),
            xgridvisible = false,
            ygridvisible = false,
        )

        _plot_dataset_records!(
            ax,
            group_records;
            errorbars=errorbars,
            xerr_col=xerr_col,
            yerr_col=yerr_col,
            logscale=logscale,
            markersize=markersize,
            strokewidth=strokewidth,
            font_scale=font_scale,
            limits=isnothing(plot_xlimits) || isnothing(limits) ? limits : plot_xlimits,
        )

        if !isnothing(limits)
            xlo, xhi = isnothing(plot_xlimits) ? limits : plot_xlimits
            ylo, yhi = limits
            limits!(ax, xlo, xhi, ylo, yhi)
        end

        icon = group_records[1].icon
        if show_icons && !isnothing(icon) && isfile(icon)
            _add_icon!(fig[row, col], icon; group_records[1].icon_kw...)
        end

        push!(axes, ax)
    end

    if length(axes) > 1
        for ax in axes[2:end]
            linkxaxes!(axes[1], ax)
            linkyaxes!(axes[1], ax)
        end
    end

    for (i, ax) in enumerate(axes)
        row, col = positions[i]
        row < nrows && hidexdecorations!(ax; grid=false)
        col > 1 && hideydecorations!(ax; grid=false)
    end

    rowgap!(fig.layout, 6)
    colgap!(fig.layout, -40)

    return fig
end

function plot(;
    RESULTDIR = Meris.DATADIR * "macro/tl-prediction/",
    category = "all",
    FILENAME = nothing,
    ext = "pdf",
    savefig::Bool = true,
    figname = nothing,
    datasets = _default_datasets(; RESULTDIR=RESULTDIR),
    min_components::Int = 40,
    kwargs...
)
    selected = if isnothing(FILENAME)
        _select_datasets(category, datasets)
    else
        [(; key=category, title=String(category), prediction_file=joinpath(RESULTDIR, FILENAME), palette=_shades(colorant"#1f77b4", 10), icon=nothing, icon_kw=(;))]
    end

    records = _plot_records(selected; min_components=min_components)
    fig = plot_C_est_vs_C_fit(records; kwargs...)

    if savefig
        suffix = isnothing(FILENAME) ? String(category) : replace(FILENAME, r"\.jld2$" => "")
        outfile = isnothing(figname) ? (Meris.FIGDIR * "tl-prediction-$(suffix).$ext") : figname
        save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

end # module TLPredictionPlotter
#/ End module
