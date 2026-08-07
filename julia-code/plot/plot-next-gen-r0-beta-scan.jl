#= Plot P(R0 > 1) against mean(beta) for selected biomass sigma_B values. =#
module NextGenR0BetaScanPlot

using CairoMakie
using DataFrames
using JLD2
using LaTeXStrings
using MakiePublication
using Printf

using Meris

const DEFAULT_INPUT = joinpath(Meris.DATADIR, "next-gen", "otu-gut1-r0-beta-scan.jld2")
const DEFAULT_OUTDIR = joinpath(Meris.FIGDIR, "next-gen")

function parse_args(args)
    options = Dict(
        "input" => DEFAULT_INPUT,
        "outdir" => DEFAULT_OUTDIR,
        "basename" => "otu-gut1-r0-beta-scan",
        "biomass-sigmas" => "",
        "beta-mean-min" => "",
        "beta-mean-max" => "",
    )

    for arg in args
        if arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code julia-code/plot/plot-next-gen-r0-beta-scan.jl [options]

            Options:
              --input=PATH             Input JLD2 scan file. Default: $(DEFAULT_INPUT)
              --outdir=PATH            Output directory. Default: $(DEFAULT_OUTDIR)
              --basename=NAME          Output filename stem. Default: otu-gut1-r0-beta-scan
              --biomass-sigmas=a,b,c   Optional comma-separated sigma_B curves to plot.
              --beta-mean-min=M        Keep only mean(beta) >= M.
              --beta-mean-max=M        Keep only mean(beta) <= M.
            """)
            exit(0)
        elseif startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "="; limit=2)
            key == "biomass-sigma" && (key = "biomass-sigmas")
            haskey(options, key) || error("Unknown option: --$key")
            options[key] = value
        else
            error("Unknown argument: $arg")
        end
    end

    parse_optional_float(value) = isempty(value) ? nothing : parse(Float64, value)
    biomass_sigmas = isempty(options["biomass-sigmas"]) ?
        nothing :
        parse.(Float64, split(options["biomass-sigmas"], ","))

    return (
        input = options["input"],
        outdir = options["outdir"],
        basename = options["basename"],
        biomass_sigmas = biomass_sigmas,
        beta_mean_min = parse_optional_float(options["beta-mean-min"]),
        beta_mean_max = parse_optional_float(options["beta-mean-max"]),
    )
end

function _load_scan(path::AbstractString)
    isfile(path) || error("Input file not found: $path")
    data = JLD2.load(path)
    haskey(data, "result") || error("Expected key `result` in $path")
    return data["result"], get(data, "parameters", nothing)
end

function biomass_group_column(result::DataFrame)
    :biomass_sigma in propertynames(result) && return :biomass_sigma
    :biomass_var in propertynames(result) && return :biomass_var
    error("Expected either `biomass_sigma` or `biomass_var` in result")
end

function _filter_result(result::DataFrame; biomass_sigmas=nothing, beta_mean_min=nothing, beta_mean_max=nothing)
    filtered = result
    group_col = biomass_group_column(filtered)
    if biomass_sigmas !== nothing
        keep = [any(isapprox(row[group_col], v; rtol=1e-8) for v in biomass_sigmas) for row in eachrow(filtered)]
        filtered = filtered[keep, :]
    end
    if beta_mean_min !== nothing
        filtered = filtered[filtered.beta_mean .>= beta_mean_min, :]
    end
    if beta_mean_max !== nothing
        filtered = filtered[filtered.beta_mean .<= beta_mean_max, :]
    end
    nrow(filtered) > 0 || error("No scan points remain after filtering")
    return filtered
end

function pow10_label(x::Real)
    x > 0 || return @sprintf("%.2g", x)

    exponent = floor(Int, log10(x))
    mantissa = x / 10.0^exponent
    if isapprox(mantissa, 1.0; rtol=1e-8, atol=1e-12)
        return "10^{$exponent}"
    end

    mantissa_text = replace(@sprintf("%.2g", mantissa), r"\.0$" => "")
    return "$(mantissa_text)\\times 10^{$exponent}"
end

function marker_indices(n::Int; max_markers::Int=18)
    n <= max_markers && return collect(1:n)
    return unique(round.(Int, range(1, n, length=max_markers)))
end

function threshold_crossing_x(x::AbstractVector, y::AbstractVector, epsilon::Real)
    idx = findfirst(>(epsilon), y)
    idx === nothing && return nothing
    idx == firstindex(y) && return x[idx]

    x0, x1 = x[idx - 1], x[idx]
    y0, y1 = y[idx - 1], y[idx]
    y1 == y0 && return x1

    t = (epsilon - y0) / (y1 - y0)
    return exp(log(x0) + t * (log(x1) - log(x0)))
end

function monotone_slopes(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n == 2 && return fill((y[2] - y[1]) / (x[2] - x[1]), 2)

    h = diff(x)
    delta = diff(y) ./ h
    m = zeros(Float64, n)
    m[1] = delta[1]
    m[end] = delta[end]

    for i in 2:(n - 1)
        if delta[i - 1] * delta[i] <= 0.0
            m[i] = 0.0
        else
            w1 = 2h[i] + h[i - 1]
            w2 = h[i] + 2h[i - 1]
            m[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i])
        end
    end

    return m
end

function smooth_marker_curve(x::AbstractVector, y::AbstractVector; points_per_interval::Int=16)
    n = length(x)
    n <= 2 && return x, y

    lx = log.(x)
    m = monotone_slopes(lx, y)
    xs = Float64[]
    ys = Float64[]

    for i in 1:(n - 1)
        h = lx[i + 1] - lx[i]
        ts = range(0.0, 1.0; length=points_per_interval + 1)
        i > 1 && (ts = Iterators.drop(ts, 1))

        for t in ts
            h00 = 2t^3 - 3t^2 + 1
            h10 = t^3 - 2t^2 + t
            h01 = -2t^3 + 3t^2
            h11 = t^3 - t^2
            push!(xs, exp(lx[i] + t * h))
            push!(ys, clamp(h00 * y[i] + h10 * h * m[i] + h01 * y[i + 1] + h11 * h * m[i + 1], 0.0, 1.0))
        end
    end

    return xs, ys
end

function plot!(parent;
    input::AbstractString=DEFAULT_INPUT,
    biomass_sigmas=nothing,
    beta_mean_min=nothing,
    beta_mean_max=nothing,
    show_legend::Bool=true,
    font_scale::Real=1.0,
    labelsize=nothing,
    ticklabelsize=nothing,
    textsize=nothing,
    markersize::Real=7,
    linewidth_scale::Real=1.0,
    legend_position=:rb,
    legend_margin=(0, 8, 2, 0),
    legend_patchlabelgap::Real=5,
    axis_kwargs=(;),
)
    result, parameters = _load_scan(input)
    result = _filter_result(
        result;
        biomass_sigmas=biomass_sigmas,
        beta_mean_min=beta_mean_min,
        beta_mean_max=beta_mean_max,
    )

    __theme = MakiePublication.theme_acs(; ishollowmarkers=[true, true])
    set_theme!(__theme)

    panel = parent isa GridLayout ? parent : GridLayout(parent)
    axis_labelsize = isnothing(labelsize) ? 17 * font_scale : labelsize
    axis_ticklabelsize = isnothing(ticklabelsize) ? 13 * font_scale : ticklabelsize
    annotation_size = isnothing(textsize) ? 15 * font_scale : textsize
    ax = Axis(
        panel[1, 1];
        xscale=log10,
        xticks=(
            [1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
            ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"],
        ),
        xlabel=L"\textrm{mean infection ratio }\bar{\beta}",
        ylabel=L"P[R_0 > 1]",
        xlabelsize=axis_labelsize,
        ylabelsize=axis_labelsize,
        xticklabelsize=axis_ticklabelsize,
        yticklabelsize=axis_ticklabelsize,
        xminorgridvisible=false,
        yminorgridvisible=false,
        axis_kwargs...,
    )

    palette = [:black, "#1f77b4", "#ff7f0e"]
    markers = [:rect, :circle, :utriangle, :diamond, :cross]
    group_col = biomass_group_column(result)
    biomass_groups = sort(unique(result[!, group_col]))
    epsilon = 1e-4
    beta_c_values = Float64[]
    beta_c_label = nothing

    for (k, biomass_group) in enumerate(biomass_groups)
        sdf = sort(result[result[!, group_col] .== biomass_group, :], :beta_mean)
        color = palette[mod1(k, length(palette))]
        marker_rows = marker_indices(nrow(sdf))
        x_marker = sdf.beta_mean[marker_rows]
        y_marker = sdf.probability_gt1[marker_rows]
        beta_c = threshold_crossing_x(sdf.beta_mean, sdf.probability_gt1, epsilon)
        if beta_c !== nothing
            push!(beta_c_values, beta_c)
            k == 1 && (beta_c_label = beta_c)
            vlines!(
                ax,
                [beta_c];
                color=color,
                linewidth=0.7 * linewidth_scale,
                linestyle=:dash,
            )
        end
        x_smooth, y_smooth = smooth_marker_curve(x_marker, y_marker)
        curve_label = if group_col == :biomass_sigma
            sigma_text = @sprintf("%.2f", biomass_group)
            latexstring("\\sigma_B=$(sigma_text)")
        else
            latexstring("\\mathrm{Var}(B)=$(pow10_label(biomass_group))")
        end
        lines!(
            ax,
            x_smooth,
            y_smooth;
            color=color,
            linewidth=1.0 * linewidth_scale,
            linestyle=:dashdot,
        )
        scatter!(
            ax,
            x_marker,
            y_marker;
            color=:white,
            strokecolor=color,
            strokewidth=1.3 * linewidth_scale,
            marker=markers[mod1(k, length(markers))],
            markersize=markersize,
            label=curve_label,
        )
    end

    if beta_c_label !== nothing
        text!(
            ax,
            1.08 * beta_c_label,
            0.62;
            text=L"\beta_c",
            rotation=0.0,
            align=(:left, :center),
            offset=(4, 0),
            fontsize=annotation_size,
            color=:black,
        )
    end

    ylims!(ax, -0.02, 1.02)
    xlims!(ax, 0.75 * minimum(result.beta_mean), maximum(result.beta_mean))
    if show_legend
        axislegend(
            ax;
            position=legend_position,
            labelsize=axis_ticklabelsize,
            patchsize=(16 * font_scale, 10 * font_scale),
            rowgap=4,
            patchlabelgap=legend_patchlabelgap,
            margin=legend_margin,
            padding=(3, 3, 3, 3),
            framevisible=false,
            backgroundcolor=(:white, 0.86),
        )
    end

    return (axis=ax, result=result, parameters=parameters)
end

function plot(;
    input::AbstractString=DEFAULT_INPUT,
    outdir::AbstractString=DEFAULT_OUTDIR,
    basename::AbstractString="otu-gut1-r0-beta-scan",
    biomass_sigmas=nothing,
    beta_mean_min=nothing,
    beta_mean_max=nothing,
    savefig::Bool=true,
)
    __theme = MakiePublication.theme_acs(; ishollowmarkers=[true, true])
    set_theme!(__theme)

    width = 1.55 * 246
    height = 0.58 * width
    fig = Figure(; size=(width, height), figure_padding=(6, 14, 6, 8))

    plot!(
        fig[1, 1];
        input=input,
        biomass_sigmas=biomass_sigmas,
        beta_mean_min=beta_mean_min,
        beta_mean_max=beta_mean_max,
    )

    if savefig
        mkpath(outdir)
        pdf_file = joinpath(outdir, "$basename.pdf")
        png_file = joinpath(outdir, "$basename.png")
        CairoMakie.save(pdf_file, fig, pt_per_unit=1)
        CairoMakie.save(png_file, fig, px_per_unit=3)
        @info "Saved R0 beta scan plot" pdf=pdf_file png=png_file
    end

    return fig
end

end # module NextGenR0BetaScanPlot

if abspath(PROGRAM_FILE) == @__FILE__
    options = NextGenR0BetaScanPlot.parse_args(ARGS)
    NextGenR0BetaScanPlot.plot(;
        input=options.input,
        outdir=options.outdir,
        basename=options.basename,
        biomass_sigmas=options.biomass_sigmas,
        beta_mean_min=options.beta_mean_min,
        beta_mean_max=options.beta_mean_max,
    )
end
