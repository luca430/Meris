#= Module to plot Taylor's-law prediction summaries for Gutenberg and GTEx data =#
#/ Start module
module TLPredictionPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings
using JLD2
using DataFrames
using Colors

#/ Modules
using Meris

const MM_TO_PT = 72.0 / 25.4
const NATURE_SINGLE_WIDTH_PT = 89.0 * MM_TO_PT
const NATURE_AXIS_LABEL_PT = 7
const NATURE_TICK_PT = 6
const NATURE_TEXT_PT = 6

const GUTENBERG_BASE = colorant"#1f77b4"
const GTEX_COLOR = colorant"#d62728"

const _GUTENBERG_HSL = convert(HSL, GUTENBERG_BASE)
const GUTENBERG_COLOR = HSL(_GUTENBERG_HSL.h, _GUTENBERG_HSL.s, 0.75)

#################
### INTERNALS ###

function _load_test_summaries(path::AbstractString)
    return JLD2.load(path, "test_summaries")
end

function _dataset_name(class)
    s = String(class)
    ls = lowercase(s)

    if startswith(ls, "guten") || ls in ("en", "it")
        return "Gutenberg"
    elseif startswith(ls, "gen-") || uppercase(s) == s
        return "GTEx"
    end

    return "Unknown"
end

function _class_name(class)
    s = String(class)
    s = replace(s, r"^gutenberg-"i => "")
    s = replace(s, r"^guten-"i => "")
    s = replace(s, r"^gen-"i => "")
    return s
end

function _dataset_color(dataset::AbstractString)
    dataset == "Gutenberg" && return GUTENBERG_COLOR
    dataset == "GTEx" && return GTEX_COLOR
    return colorant"#555555"
end

function _plot_order(classes)
    return sort(
        collect(classes);
        by = c -> (_dataset_name(c) == "GTEx" ? 0 : 1, String(c))
    )
end

#################
### FUNCTIONS ###

function plot_C_est_vs_C_fit(
    test_summaries::Dict;
    errorbars::Bool = true,
    xerr_col::Symbol = :C_est_err,
    yerr_col::Symbol = :C_fit_err,
    logscale::Bool = true,
    markersize::Int = 6,
    strokewidth::Float64 = 0.7,
    font_scale::Float64 = 1.3,
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    markers = [:circle, :rect, :diamond, :utriangle, :dtriangle, :pentagon, :hexagon, :star4]
    classes = _plot_order(keys(test_summaries))

    fig = Figure(
        size = (NATURE_SINGLE_WIDTH_PT, 0.88 * NATURE_SINGLE_WIDTH_PT),
        figure_padding = (4, 10, 4, 4)
    )

    ax = Axis(
        fig[1, 1],
        xlabel = L"C_{\mathrm{est}}",
        ylabel = L"C_{\mathrm{fit}}",
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

    for (j, class) in enumerate(classes)
        df = test_summaries[class]
        d = copy(df)
        dataset = _dataset_name(class)
        color = _dataset_color(dataset)
        marker = markers[mod1(j, length(markers))]

        filter!(
            [:C_est, :C_fit, xerr_col, yerr_col] =>
                (x, y, xe, ye) ->
                    all(isfinite, (x, y, xe, ye)) &&
                    x > 0 && y > 0 &&
                    x - xe > 0 &&
                    y - ye > 0,
            d
        )

        if nrow(d) == 0
            @warn "Skipping class $class: no valid points"
            continue
        end

        x = d.C_est
        y = d.C_fit
        xerr = d[:, xerr_col]
        yerr = d[:, yerr_col]

        append!(all_x, x)
        append!(all_y, y)

        if errorbars
            errorbars!(
                ax,
                x,
                y,
                yerr;
                color = color,
            whiskerwidth = 4,
            linewidth = 0.7,
            )

            errorbars!(
                ax,
                x,
                y,
                xerr;
                direction = :x,
                color = color,
                whiskerwidth = 4,
                linewidth = 0.7,
            )
        end

        scatter!(
            ax,
            x,
            y;
            label = "$(dataset): $(_class_name(class))",
            marker = marker,
            color = (:white, 1.0),
            strokecolor = color,
            markersize = markersize,
            strokewidth = strokewidth,
        )
    end

    if !isempty(all_x) && !isempty(all_y)
        lo = minimum(vcat(all_x, all_y))
        hi = maximum(vcat(all_x, all_y))

        xs = logscale ? exp10.(range(log10(lo), log10(hi); length=300)) :
                        range(lo, hi; length=300)

        lines!(
            ax,
            xs,
            xs;
            linestyle = (:dash, :dense),
            color = :black,
            linewidth = 1.2,
            label = L"C_{\mathrm{fit}} = C_{\mathrm{est}}",
        )

        limits!(ax, lo, hi, lo, hi)
    end

    axislegend(
        ax;
        position = :lt,
        patchsize = (8, 8),
        rowgap = 0,
        labelsize = NATURE_TEXT_PT * font_scale,
        padding = 2,
    )

    return fig
end

function plot(;
    RESULTDIR = Meris.DATADIR * "macro/tl-prediction/",
    category = "linguistic",
    FILENAME = nothing,
    ext = "pdf",
    savefig::Bool = true,
    figname = nothing,
    kwargs...
)
    file = isnothing(FILENAME) ? "$(category).jld2" : FILENAME
    test_summaries = _load_test_summaries(RESULTDIR * file)
    fig = plot_C_est_vs_C_fit(test_summaries; kwargs...)

    if savefig
        outfile = isnothing(figname) ? (Meris.FIGDIR * "tl-prediction-$(category).$ext") : figname
        save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

end # module TLPredictionPlotter
#/ End module
