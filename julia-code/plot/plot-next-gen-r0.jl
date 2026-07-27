#= Plot the sample-level R0 distribution from next-generation matrix analysis. =#
module NextGenR0Plot

using CairoMakie
using DataFrames
using JLD2
using LaTeXStrings
using MakiePublication
using Printf
using Statistics

using Meris

const DEFAULT_INPUT = joinpath(Meris.DATADIR, "next-gen", "otu-gut1-r0.jld2")
const DEFAULT_OUTDIR = joinpath(Meris.FIGDIR, "next-gen")

function parse_args(args)
    options = Dict(
        "input" => DEFAULT_INPUT,
        "outdir" => DEFAULT_OUTDIR,
        "basename" => "otu-gut1-log-r0-distribution",
        "nbins" => "28",
        "distribution-ids" => "",
    )

    for arg in args
        if arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code julia-code/plot/plot-next-gen-r0.jl [options]

            Options:
              --input=PATH              Input JLD2 file. Default: $(DEFAULT_INPUT)
              --outdir=PATH             Output directory. Default: $(DEFAULT_OUTDIR)
              --basename=NAME           Output filename stem. Default: otu-gut1-log-r0-distribution
              --nbins=N                 Number of bins. Default: 28
              --distribution-ids=a,b     Optional comma-separated subset of distribution IDs.
            """)
            exit(0)
        elseif startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "="; limit=2)
            key == "distribution-id" && (key = "distribution-ids")
            haskey(options, key) || error("Unknown option: --$key")
            options[key] = value
        else
            error("Unknown argument: $arg")
        end
    end

    distribution_ids = isempty(options["distribution-ids"]) ?
        nothing :
        parse.(Int, split(options["distribution-ids"], ","))

    return (
        input = options["input"],
        outdir = options["outdir"],
        basename = options["basename"],
        nbins = parse(Int, replace(options["nbins"], "_" => "")),
        distribution_ids = distribution_ids,
    )
end

function _load_r0(path::AbstractString)
    isfile(path) || error("Input file not found: $path")
    data = JLD2.load(path)
    haskey(data, "r0df") || error("Expected key `r0df` in $path")
    return data["r0df"], get(data, "parameters", nothing), get(data, "distribution_parameters", nothing)
end

function _prepare_r0df!(r0df::DataFrame)
    if :distribution_id ∉ propertynames(r0df)
        r0df.distribution_id = fill(1, nrow(r0df))
    end
    return r0df
end

function _filter_distributions(r0df::DataFrame, distribution_ids)
    distribution_ids === nothing && return r0df
    filtered = r0df[in.(r0df.distribution_id, Ref(distribution_ids)), :]
    nrow(filtered) > 0 || error("No rows found for distribution IDs $(join(distribution_ids, ","))")
    return filtered
end

function _distribution_label(distribution_id::Int, distribution_params)
    distribution_params === nothing && return "distribution $distribution_id"
    matches = distribution_params[distribution_params.distribution_id .== distribution_id, :]
    nrow(matches) == 0 && return "distribution $distribution_id"

    row = first(eachrow(matches))
    return latexstring("\\bar{B}=$(pow10_label(row.biomass_mean)),\\;\\sigma^2_B=$(pow10_label(row.biomass_var))")
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

function _hist_density(values::AbstractVector, edges::AbstractVector)
    counts = zeros(Float64, length(edges) - 1)
    for value in values
        idx = searchsortedlast(edges, value)
        idx = clamp(idx, 1, length(counts))
        counts[idx] += 1
    end

    widths = diff(edges)
    density = counts ./ (sum(counts) .* widths)
    return density
end

function _stairs_xy(edges::AbstractVector, density::AbstractVector)
    x = Float64[]
    y = Float64[]
    for i in eachindex(density)
        push!(x, edges[i], edges[i + 1])
        push!(y, density[i], density[i])
    end
    return x, y
end

function plot(;
    input::AbstractString=DEFAULT_INPUT,
    outdir::AbstractString=DEFAULT_OUTDIR,
    basename::AbstractString="otu-gut1-log-r0-distribution",
    nbins::Int=28,
    distribution_ids=nothing,
    savefig::Bool=true,
)
    r0df, parameters, distribution_params = _load_r0(input)
    _prepare_r0df!(r0df)
    r0df = _filter_distributions(r0df, distribution_ids)
    r0df = r0df[r0df.R0 .> 0.0, :]
    nrow(r0df) > 0 || error("No positive R0 values found in $input")
    r0df.logR0 = log10.(r0df.R0)

    __theme = MakiePublication.theme_acs(; ishollowmarkers=[true, true])
    set_theme!(__theme)

    width = 1.15 * 246
    height = 3 * width / 4.67
    fig = Figure(; size=(width, height), figure_padding=(6, 8, 4, 12))
    ax = Axis(
        fig[1, 1],
        xlabel=L"\log_{10} R_0",
        ylabel=L"\textrm{density}",
        xlabelsize=11,
        ylabelsize=11,
    )

    log_values = r0df.logR0
    lo, hi = extrema(log_values)
    lo == hi && ((lo, hi) = (lo - 0.5, hi + 0.5))
    edges = collect(range(lo, hi, length=nbins + 1))
    dist_ids = sort(unique(r0df.distribution_id))
    palette = Makie.wong_colors()
    xlims!(ax, lo, hi)

    if hi > 0.0
        vspan!(
            ax,
            0.0,
            hi;
            color=(:lightcoral, 0.16),
        )
    end

    for (k, distribution_id) in enumerate(dist_ids)
        sdf = r0df[r0df.distribution_id .== distribution_id, :]
        density = _hist_density(sdf.logR0, edges)
        x, y = _stairs_xy(edges, density)
        color = palette[mod1(k, length(palette))]

        lines!(
            ax,
            x,
            y;
            color=color,
            linewidth=1.5,
            label=_distribution_label(distribution_id, distribution_params),
        )
        vlines!(ax, [mean(sdf.logR0)]; color=color, linewidth=0.8, linestyle=:dash)
    end

    axislegend(
        ax;
        position=:lt,
        labelsize=8,
        patchsize=(10, 6),
        rowgap=0,
        padding=(2, 2, 2, 2),
    )

    title = "OTU $(first(r0df.class))"
    if parameters !== nothing
        if :n_distributions in propertynames(parameters)
            title *= ", distributions=$(parameters.n_distributions), n=$(parameters.ncomponents)"
        elseif :connectivity in propertynames(parameters)
            title *= ", c=$(parameters.connectivity), n=$(parameters.ncomponents)"
        end
    end
    Label(fig[0, 1], title; fontsize=10, tellwidth=false)

    if savefig
        mkpath(outdir)
        pdf_file = joinpath(outdir, "$basename.pdf")
        png_file = joinpath(outdir, "$basename.png")
        CairoMakie.save(pdf_file, fig, pt_per_unit=1)
        CairoMakie.save(png_file, fig, px_per_unit=3)
        @info "Saved R0 distribution plot" pdf=pdf_file png=png_file
    end

    return fig
end

end # module NextGenR0Plot

if abspath(PROGRAM_FILE) == @__FILE__
    options = NextGenR0Plot.parse_args(ARGS)
    NextGenR0Plot.plot(;
        input=options.input,
        outdir=options.outdir,
        basename=options.basename,
        nbins=options.nbins,
        distribution_ids=options.distribution_ids,
    )
end
