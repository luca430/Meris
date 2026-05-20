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

#/ Modules
using Meris

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
        vcat(_shades(bases[4], 10)[1:7], _shades(bases[5], 8)),
    ]

    return [
        (;
            key="linguistic",
            title="Linguistic",
            prediction_file=joinpath(RESULTDIR, "linguistic.jld2"),
            palette=palettes[1],
        ),
        (;
            key="microbial",
            title="Microbial",
            prediction_file=joinpath(RESULTDIR, "microbial.jld2"),
            palette=palettes[2],
        ),
        (;
            key="social",
            title="Social",
            prediction_file=joinpath(RESULTDIR, "social.jld2"),
            palette=palettes[3],
        ),
        (;
            key="biology",
            title="Biology",
            prediction_file=joinpath(RESULTDIR, "biology.jld2"),
            palette=palettes[4],
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

function _plot_records(datasets; min_components::Int=11)
    records = NamedTuple[]

    for spec in datasets
        summaries, component_bins, summary_kind = _load_prediction_data(spec.prediction_file)
        classes = sort(collect(keys(summaries)); by=string)

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

#################
### FUNCTIONS ###

function plot_C_est_vs_C_fit(
    records;
    errorbars::Bool = false,
    xerr_col::Symbol = :C_est_err,
    yerr_col::Symbol = :C_fit_err,
    logscale::Bool = true,
    markersize::Int = 6,
    strokewidth::Float64 = 0.7,
    font_scale::Float64 = 2.0,
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    fig = Figure(
        size = (300, 256),
        figure_padding = (4, 6, 4, 4)
    )

    ax = Axis(
        fig[1, 1],
        xlabel = L"\Omega_{\mathrm{fit},A}",
        ylabel = L"\Omega_{\mathrm{fit},B}",
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

    all_x = Float64[]
    all_y = Float64[]
    fit_groups = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()

    for record in records
        d = copy(record.df)
        color = record.color

        filter!(
            [:C_est, :C_fit] =>
                (x, y) -> all(isfinite, (x, y)) && x > 0 && y > 0,
            d
        )

        if nrow(d) == 0
            @warn "Skipping class $(record.class): no valid points"
            continue
        end

        x = d.C_est
        y = d.C_fit

        append!(all_x, x)
        append!(all_y, y)

        fit_key = "$(record.dataset_key):$(record.class)"
        gx, gy = get!(fit_groups, fit_key, (Float64[], Float64[]))
        append!(gx, x)
        append!(gy, y)

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
                    color = color,
                    whiskerwidth = 4,
                    linewidth = 0.7,
                )

                errorbars!(
                    ax,
                    d_err.C_est,
                    d_err.C_fit,
                    d_err[:, xerr_col];
                    direction = :x,
                    color = color,
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
            strokecolor = color,
            markersize = markersize,
            strokewidth = strokewidth,
        )
    end

    if !isempty(all_x) && !isempty(all_y)
        slopes = [
            sum(x .* y) / sum(abs2, x)
            for (x, y) in values(fit_groups)
            if !isempty(x) && sum(abs2, x) > 0
        ]
        a_fit = mean(slopes)
        a_fit_err = length(slopes) > 1 ? std(slopes) / sqrt(length(slopes)) : 0.0
        fit_y = a_fit .* all_x
        lo = minimum(vcat(all_x, all_y, fit_y))
        hi = maximum(vcat(all_x, all_y, fit_y))

        xs = logscale ? exp10.(range(log10(lo), log10(hi); length=300)) :
                        range(lo, hi; length=300)

        lines!(
            ax,
            xs,
            a_fit .* xs;
            linestyle = (:dash, :dense),
            color = :black,
            linewidth = 2.0,
            label = L"\Omega_{\mathrm{fit},B} = a\Omega_{\mathrm{fit},A}",
        )

        text!(
            ax,
            0.05,
            0.89;
            text = L"a = %$(round(a_fit; digits=2)) \pm %$(round(a_fit_err; digits=2))",
            space = :relative,
            align = (:left, :top),
            fontsize = NATURE_TEXT_PT * font_scale,
            color = :black,
        )

        limits!(ax, lo, hi, lo, hi)
    end

    axislegend(
        ax;
        position = :lt,
        patchsize = (14, 8),
        nbanks = 2,
        labelsize = NATURE_TEXT_PT * font_scale,
        padding = 2,
        framevisible = false,
    )

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
    min_components::Int = 50,
    kwargs...
)
    selected = if isnothing(FILENAME)
        _select_datasets(category, datasets)
    else
        [(; key=category, title=String(category), prediction_file=joinpath(RESULTDIR, FILENAME), palette=_shades(colorant"#1f77b4", 10))]
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
