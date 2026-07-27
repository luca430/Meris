#= Plot P(R0 > 1) heatmap from next-generation matrix grid analysis. =#
module NextGenR0GridPlot

using CairoMakie
using DataFrames
using JLD2
using LaTeXStrings
using MakiePublication

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
        "biomass-var-min" => "",
        "biomass-var-max" => "",
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
              --log-color=true|false   Plot log10 probabilities; zero cells are grey. Default: false
              --biomass-var-min=V      Keep only Var(B) >= V.
              --biomass-var-max=V      Keep only Var(B) <= V.
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
        biomass_var_min = parse_optional_float(options["biomass-var-min"]),
        biomass_var_max = parse_optional_float(options["biomass-var-max"]),
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

function _filter_result(result::DataFrame; biomass_var_min=nothing, biomass_var_max=nothing,
                        beta_mean_min=nothing, beta_mean_max=nothing)
    filtered = result
    if biomass_var_min !== nothing
        filtered = filtered[filtered.biomass_var .>= biomass_var_min, :]
    end
    if biomass_var_max !== nothing
        filtered = filtered[filtered.biomass_var .<= biomass_var_max, :]
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

function _grid_matrix(result::DataFrame, biomass_vars, beta_means)
    biomass_vars = sort(unique(result.biomass_var))
    beta_means = sort(unique(result.beta_mean))

    z = fill(NaN, length(biomass_vars), length(beta_means))
    lookup = Dict(
        (row.biomass_var, row.beta_mean) => row.probability_gt1
        for row in eachrow(result)
    )

    for i in eachindex(biomass_vars), j in eachindex(beta_means)
        z[i, j] = lookup[(biomass_vars[i], beta_means[j])]
    end

    return collect(biomass_vars), collect(beta_means), z
end

function _color_values(probabilities; log_color::Bool=false, pmin=nothing, pmax=nothing)
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

    colorrange = (
        isnothing(pmin) ? 0.0 : pmin,
        isnothing(pmax) ? 1.0 : pmax,
    )
    return probabilities, colorrange, L"P(R_0 > 1)"
end

function plot(;
    input::AbstractString=DEFAULT_INPUT,
    outdir::AbstractString=DEFAULT_OUTDIR,
    basename::AbstractString="otu-gut1-r0-grid-heatmap",
    pmin=nothing,
    pmax=nothing,
    log_color::Bool=false,
    biomass_var_min=nothing,
    biomass_var_max=nothing,
    beta_mean_min=nothing,
    beta_mean_max=nothing,
    savefig::Bool=true,
)
    result, biomass_vars, beta_means, parameters = _load_grid(input)
    result = _filter_result(
        result;
        biomass_var_min=biomass_var_min,
        biomass_var_max=biomass_var_max,
        beta_mean_min=beta_mean_min,
        beta_mean_max=beta_mean_max,
    )
    biomass_vars, beta_means, probabilities = _grid_matrix(result, biomass_vars, beta_means)
    color_values, colorrange, colorbar_label = _color_values(
        probabilities;
        log_color=log_color,
        pmin=pmin,
        pmax=pmax,
    )

    __theme = MakiePublication.theme_acs(; ishollowmarkers=[true, true])
    set_theme!(__theme)

    width = 1.2 * 246
    height = 0.95 * width
    fig = Figure(; size=(width, height), figure_padding=(8, 12, 6, 12))

    ax = Axis(
        fig[1, 1],
        xlabel=L"\log_{10}\,\mathrm{Var}(B)",
        ylabel=L"\log_{10}\,\langle \beta \rangle",
        xlabelsize=11,
        ylabelsize=11,
    )

    hm = heatmap!(
        ax,
        log10.(biomass_vars),
        log10.(beta_means),
        color_values;
        colormap=:viridis,
        colorrange=colorrange,
        nan_color=:lightgray,
    )

    Colorbar(
        fig[1, 2],
        hm;
        label=colorbar_label,
        labelsize=10,
        ticklabelsize=8,
        width=10,
    )

    title = "OTU GUT1"
    if parameters !== nothing
        title *= ", c=$(parameters.connectivity), runs=$(parameters.n_runs)"
    end
    log_color && (title *= ", log color")
    Label(fig[0, 1:2], title; fontsize=10, tellwidth=false)

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
        biomass_var_min=options.biomass_var_min,
        biomass_var_max=options.biomass_var_max,
        beta_mean_min=options.beta_mean_min,
        beta_mean_max=options.beta_mean_max,
    )
end
