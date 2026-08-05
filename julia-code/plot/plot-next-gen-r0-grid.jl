#= Plot P(R0 > 1) heatmap from next-generation matrix grid analysis. =#
module NextGenR0GridPlot

using CairoMakie
using DataFrames
using JLD2
using LaTeXStrings
using MakiePublication
using Makie

using Meris

const DEFAULT_INPUT = joinpath(Meris.DATADIR, "next-gen", "otu-gut1-r0-grid.jld2")
const DEFAULT_OUTDIR = joinpath(Meris.FIGDIR, "next-gen")

function parse_args(args)
    options = Dict(
        "input" => DEFAULT_INPUT,
        "outdir" => DEFAULT_OUTDIR,
        "basename" => "otu-gut1-r0-grid-heatmap",
        "pmin" => "",
        "pmax" => "",
        "log-color" => "false",
        "probability-gamma" => "0.35",
        "biomass-sigma-min" => "",
        "biomass-sigma-max" => "",
        "beta-mean-min" => "",
        "beta-mean-max" => "",
    )

    for arg in args
        if arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code julia-code/plot/plot-next-gen-r0-grid.jl [options]

            Options:
              --input=PATH             Input JLD2 grid file. Default: $(DEFAULT_INPUT)
              --outdir=PATH            Output directory. Default: $(DEFAULT_OUTDIR)
              --basename=NAME          Output filename stem. Default: otu-gut1-r0-grid-heatmap
              --pmin=P                 Heatmap color minimum. Default: 0, or auto for --log-color=true
              --pmax=P                 Heatmap color maximum. Default: 1, or auto for --log-color=true
              --log-color=true|false   Plot log10 probabilities; zero cells are white. Default: false
              --probability-gamma=G    Display P^G for linear color; smaller emphasizes low P. Default: 0.35
              --biomass-sigma-min=S    Keep only sigma_B >= S.
              --biomass-sigma-max=S    Keep only sigma_B <= S.
              --beta-mean-min=M        Keep only mean(beta) >= M.
              --beta-mean-max=M        Keep only mean(beta) <= M.
            """)
            exit(0)
        elseif startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "="; limit=2)
            haskey(options, key) || error("Unknown option: --$key")
            options[key] = value
        else
            error("Unknown argument: $arg")
        end
    end

    parse_optional_float(value) = isempty(value) ? nothing : parse(Float64, value)
    return (
        input = options["input"],
        outdir = options["outdir"],
        basename = options["basename"],
        pmin = parse_optional_float(options["pmin"]),
        pmax = parse_optional_float(options["pmax"]),
        log_color = lowercase(options["log-color"]) in ("true", "1", "yes"),
        probability_gamma = parse(Float64, options["probability-gamma"]),
        biomass_sigma_min = parse_optional_float(options["biomass-sigma-min"]),
        biomass_sigma_max = parse_optional_float(options["biomass-sigma-max"]),
        beta_mean_min = parse_optional_float(options["beta-mean-min"]),
        beta_mean_max = parse_optional_float(options["beta-mean-max"]),
    )
end

function _load_grid(path::AbstractString)
    isfile(path) || error("Input file not found: $path")
    data = JLD2.load(path)
    haskey(data, "result") || error("Expected key `result` in $path")
    return data["result"], get(data, "biomass_vars", nothing), get(data, "beta_means", nothing), get(data, "parameters", nothing)
end

function biomass_axis_column(result::DataFrame)
    :biomass_sigma in propertynames(result) && return :biomass_sigma
    :biomass_var in propertynames(result) && return :biomass_var
    error("Expected either `biomass_sigma` or `biomass_var` in result")
end

function _filter_result(result::DataFrame; biomass_sigma_min=nothing, biomass_sigma_max=nothing,
                        beta_mean_min=nothing, beta_mean_max=nothing)
    filtered = result
    axis_col = biomass_axis_column(filtered)
    if biomass_sigma_min !== nothing
        filtered = filtered[filtered[!, axis_col] .>= biomass_sigma_min, :]
    end
    if biomass_sigma_max !== nothing
        filtered = filtered[filtered[!, axis_col] .<= biomass_sigma_max, :]
    end
    if beta_mean_min !== nothing
        filtered = filtered[filtered.beta_mean .>= beta_mean_min, :]
    end
    if beta_mean_max !== nothing
        filtered = filtered[filtered.beta_mean .<= beta_mean_max, :]
    end
    nrow(filtered) > 0 || error("No grid points remain after filtering")
    return filtered
end

function _grid_matrix(result::DataFrame, biomass_values, beta_means)
    axis_col = biomass_axis_column(result)
    biomass_values = sort(unique(result[!, axis_col]))
    beta_means = sort(unique(result.beta_mean))

    z = fill(NaN, length(biomass_values), length(beta_means))
    lookup = Dict(
        (row[axis_col], row.beta_mean) => row.probability_gt1
        for row in eachrow(result)
    )

    for i in eachindex(biomass_values), j in eachindex(beta_means)
        z[i, j] = lookup[(biomass_values[i], beta_means[j])]
    end

    return collect(biomass_values), collect(beta_means), z, axis_col
end

function _color_values(probabilities; log_color::Bool=false, pmin=nothing, pmax=nothing,
                       probability_gamma::Real=1.0)
    if log_color
        values = copy(probabilities)
        values[values .<= 0.0] .= NaN
        values = log10.(values)
        finite_values = values[isfinite.(values)]
        isempty(finite_values) && error("No positive probabilities available for log-color plot")
        colorrange = (
            isnothing(pmin) ? minimum(finite_values) : pmin,
            isnothing(pmax) ? maximum(finite_values) : pmax,
        )
        return values, colorrange, L"\log_{10} P(R_0 > 1)"
    end

    probability_gamma > 0.0 || error("--probability-gamma must be positive")
    colorrange = (
        isnothing(pmin) ? 0.0 : pmin^probability_gamma,
        isnothing(pmax) ? 1.0 : pmax^probability_gamma,
    )
    return probabilities .^ probability_gamma, colorrange, L"P(R_0 > 1)"
end

function _colorbar_ticks(; log_color::Bool=false, probability_gamma::Real=1.0)
    log_color && return Makie.automatic
    probabilities = [0.0, 1e-2, 0.1, 0.5, 1.0]
    labels = [L"0", L"10^{-2}", L"0.1", L"0.5", L"1.0"]
    return (probabilities .^ probability_gamma, labels)
end

function plot(;
    input::AbstractString=DEFAULT_INPUT,
    outdir::AbstractString=DEFAULT_OUTDIR,
    basename::AbstractString="otu-gut1-r0-grid-heatmap",
    pmin=nothing,
    pmax=nothing,
    log_color::Bool=false,
    probability_gamma::Real=0.35,
    biomass_sigma_min=nothing,
    biomass_sigma_max=nothing,
    beta_mean_min=nothing,
    beta_mean_max=nothing,
    savefig::Bool=true,
)
    result, biomass_vars, beta_means, parameters = _load_grid(input)
    result = _filter_result(
        result;
        biomass_sigma_min=biomass_sigma_min,
        biomass_sigma_max=biomass_sigma_max,
        beta_mean_min=beta_mean_min,
        beta_mean_max=beta_mean_max,
    )
    biomass_values, beta_means, probabilities, biomass_axis = _grid_matrix(result, biomass_vars, beta_means)
    color_values, colorrange, colorbar_label = _color_values(
        probabilities;
        log_color=log_color,
        pmin=pmin,
        pmax=pmax,
        probability_gamma=probability_gamma,
    )

    __theme = MakiePublication.theme_acs(; ishollowmarkers=[true, true])
    set_theme!(__theme)

    width = 1.2 * 246
    height = 0.95 * width
    fig = Figure(; size=(width, height), figure_padding=(6, 14, 6, 8))

    ax = Axis(
        fig[1, 1],
        xlabel=biomass_axis == :biomass_sigma ? L"\sigma_B" : L"\log_{10}\,\mathrm{Var}(B)",
        ylabel=L"\log_{10}\,\langle \beta \rangle",
        xlabelsize=17,
        ylabelsize=17,
        xticklabelsize=13,
        yticklabelsize=13,
        xminorgridvisible=false,
        yminorgridvisible=false,
    )

    blue_colormap = cgrad([:white, "#d8edf8", "#75b4d8", "#1f77b4"])
    hm = heatmap!(
        ax,
        biomass_axis == :biomass_sigma ? biomass_values : log10.(biomass_values),
        log10.(beta_means),
        color_values;
        colormap=blue_colormap,
        colorrange=colorrange,
        interpolate=true,
        nan_color=:white,
    )

    if !log_color && minimum(probabilities) <= 1e-3 <= maximum(probabilities)
        contour!(
            ax,
            biomass_axis == :biomass_sigma ? biomass_values : log10.(biomass_values),
            log10.(beta_means),
            probabilities;
            levels=[1e-3],
            color=:black,
            linewidth=1.2,
            linestyle=:dash,
        )
        xlo, xhi = extrema(biomass_axis == :biomass_sigma ? biomass_values : log10.(biomass_values))
        ylo, yhi = extrema(log10.(beta_means))
        text!(
            ax,
            xlo + 0.50 * (xhi - xlo),
            ylo + 0.21 * (yhi - ylo);
            text=L"P \approx 0",
            color=:black,
            fontsize=15,
            rotation=0.10pi,
            align=(:center, :center),
        )
    end

    if !log_color && minimum(probabilities) <= 0.95 <= maximum(probabilities)
        contour!(
            ax,
            biomass_axis == :biomass_sigma ? biomass_values : log10.(biomass_values),
            log10.(beta_means),
            probabilities;
            levels=[0.95],
            color=:white,
            linewidth=1.4,
            linestyle=:dash,
        )
        xlo, xhi = extrema(biomass_axis == :biomass_sigma ? biomass_values : log10.(biomass_values))
        ylo, yhi = extrema(log10.(beta_means))
        text!(
            ax,
            xlo + 0.54 * (xhi - xlo),
            ylo + 0.83 * (yhi - ylo);
            text=L"P \approx 1",
            color=:white,
            fontsize=16,
            rotation=0.18pi,
            align=(:center, :center),
        )
    end

    Colorbar(
        fig[1, 2],
        hm;
        label=colorbar_label,
        ticks=_colorbar_ticks(; log_color=log_color, probability_gamma=probability_gamma),
        labelsize=14,
        ticklabelsize=12,
        width=8,
    )

    if savefig
        mkpath(outdir)
        pdf_file = joinpath(outdir, "$basename.pdf")
        png_file = joinpath(outdir, "$basename.png")
        CairoMakie.save(pdf_file, fig, pt_per_unit=1)
        CairoMakie.save(png_file, fig, px_per_unit=3)
        @info "Saved R0 grid heatmap" pdf=pdf_file png=png_file
    end

    return fig
end

end # module NextGenR0GridPlot

if abspath(PROGRAM_FILE) == @__FILE__
    options = NextGenR0GridPlot.parse_args(ARGS)
    NextGenR0GridPlot.plot(;
        input=options.input,
        outdir=options.outdir,
        basename=options.basename,
        pmin=options.pmin,
        pmax=options.pmax,
        log_color=options.log_color,
        probability_gamma=options.probability_gamma,
        biomass_sigma_min=options.biomass_sigma_min,
        biomass_sigma_max=options.biomass_sigma_max,
        beta_mean_min=options.beta_mean_min,
        beta_mean_max=options.beta_mean_max,
    )
end
