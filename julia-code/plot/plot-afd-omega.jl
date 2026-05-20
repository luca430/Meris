#= Plot downsampled AFD panels for components in a selected Omega bin. =#

include(joinpath(@__DIR__, "plot-afd.jl"))

using DataFrames
using Meris
using .AFDPlotter

const DEFAULT_OCCUPANCY = 0.1
const DEFAULT_MEAN_OMEGA_THRESHOLD = 10.0
const DEFAULT_LOG10_HALF_WIDTH = 0.25
const DEFAULT_RESCALE_BY_OCCUPANCY = false

function _omega_label_value(omega::Real)
    logomega = log10(Float64(omega))
    if isapprox(logomega, round(logomega); atol=1e-8)
        return "1e$(Int(round(logomega)))"
    end
    return replace(string(round(Float64(omega), sigdigits=4)), "." => "p")
end

function _component_omega_stats(df::DataFrame)
    nsamples = length(unique(df.sample_id))
    stats = DataFrame(
        component_id=eltype(df.component_id)[],
        occupancy=Float64[],
        mean=Float64[],
        variance=Float64[],
        omega=Float64[],
        mean_omega=Float64[],
    )

    for sdf in groupby(df, :component_id)
        counts = Float64.(sdf.counts)
        μ = sum(counts) / nsamples
        σ2 = sum(abs2, counts) / nsamples - μ^2

        (!isfinite(μ) || !isfinite(σ2) || μ <= 0) && continue

        Ω = (σ2 - μ) / μ^2
        (!isfinite(Ω) || Ω <= 0) && continue

        push!(stats, (;
            component_id=sdf.component_id[1],
            occupancy=nrow(sdf) / nsamples,
            mean=μ,
            variance=σ2,
            omega=Ω,
            mean_omega=μ * Ω,
        ))
    end

    return stats
end

function _omega_value_filter(
    df::DataFrame;
    occ::Float64=DEFAULT_OCCUPANCY,
    omega::Real,
    log10_half_width::Float64=DEFAULT_LOG10_HALF_WIDTH,
    mean_omega_threshold::Float64=DEFAULT_MEAN_OMEGA_THRESHOLD,
)
    nrow(df) == 0 && return df
    omega > 0 || throw(ArgumentError("omega must be positive"))
    log10_half_width > 0 || throw(ArgumentError("log10_half_width must be positive"))

    target = log10(Float64(omega))
    stats = _component_omega_stats(df)
    filter!(
        row ->
            row.occupancy > occ &&
            row.mean_omega > mean_omega_threshold &&
            abs(log10(row.omega) - target) <= log10_half_width,
        stats,
    )

    keep = Set(stats.component_id)
    return filter(row -> row.component_id in keep, df)
end

function plot_downsampled_omega_afd(;
    omega::Real,
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    occ::Float64=DEFAULT_OCCUPANCY,
    mean_omega_threshold::Float64=DEFAULT_MEAN_OMEGA_THRESHOLD,
    log10_half_width::Float64=DEFAULT_LOG10_HALF_WIDTH,
    rescale_by_occupancy::Bool=DEFAULT_RESCALE_BY_OCCUPANCY,
    downsampled_dir=joinpath(Meris.DATADIR, "downsampled"),
    kwargs...
)
    omega > 0 || throw(ArgumentError("omega must be positive"))
    outfile = if isnothing(figname)
        Meris.FIGDIR * "afd-downsampled-omega-$(_omega_label_value(omega)).$ext"
    else
        figname
    end

    return AFDPlotter.plot_afd(;
        ext=ext,
        savefig=savefig,
        figname=outfile,
        datasets=AFDPlotter._default_downsampled_datasets(; downsampled_dir=downsampled_dir),
        component_filter=(df; occ) -> _omega_value_filter(
            df;
            occ=occ,
            omega=omega,
            log10_half_width=log10_half_width,
            mean_omega_threshold=mean_omega_threshold,
        ),
        occ=occ,
        rescale_by_occupancy=rescale_by_occupancy,
        kwargs...
    )
end

function _usage()
    return """
    Usage:
      julia --project=julia-code julia-code/plot/plot-afd-omega.jl --omega=VALUE [options]
      julia --project=julia-code julia-code/plot/plot-afd-omega.jl --omega-exponent=K [options]

    Options:
      --omega=VALUE                Target Omega value.
      --omega-exponent=K           Target Omega value as 10^K.
      --log10-half-width=VALUE     Half-width of the log10(Omega) bin. Default: $(DEFAULT_LOG10_HALF_WIDTH)
      --occ=VALUE                  Occupancy threshold. Default: $(DEFAULT_OCCUPANCY)
      --mean-omega-threshold=VALUE Mean_i * Omega_i threshold. Default: $(DEFAULT_MEAN_OMEGA_THRESHOLD)
      --rescale-by-occupancy=BOOL  Rescale AFD mean/variance by occupancy. Default: $(DEFAULT_RESCALE_BY_OCCUPANCY)
      --ext=EXT                    Output extension. Default: pdf
      --figname=PATH               Output path. Default: Meris.FIGDIR/afd-downsampled-omega-<omega>.EXT
      --downsampled-dir=DIR        Directory containing grouped downsampled .jld2 files.
    """
end

function _parse_bool(value::AbstractString)
    v = lowercase(strip(value))
    v in ("true", "t", "yes", "y", "1") && return true
    v in ("false", "f", "no", "n", "0") && return false
    throw(ArgumentError("Expected a boolean value, got: $value"))
end

function _parse_cli(args)
    opts = Dict{String,String}()
    for arg in args
        if arg in ("-h", "--help")
            println(_usage())
            exit(0)
        elseif startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "="; limit=2)
            opts[key] = value
        else
            error("Unknown argument: $arg\n$(_usage())")
        end
    end

    has_omega = haskey(opts, "omega")
    has_exponent = haskey(opts, "omega-exponent")
    has_omega || has_exponent || error("Missing --omega or --omega-exponent\n$(_usage())")
    !(has_omega && has_exponent) || error("Use only one of --omega and --omega-exponent")

    omega = has_omega ? parse(Float64, opts["omega"]) : 10.0 ^ parse(Float64, opts["omega-exponent"])
    return (;
        omega=omega,
        ext=get(opts, "ext", "pdf"),
        figname=get(opts, "figname", nothing),
        occ=parse(Float64, get(opts, "occ", string(DEFAULT_OCCUPANCY))),
        mean_omega_threshold=parse(Float64, get(opts, "mean-omega-threshold", string(DEFAULT_MEAN_OMEGA_THRESHOLD))),
        log10_half_width=parse(Float64, get(opts, "log10-half-width", string(DEFAULT_LOG10_HALF_WIDTH))),
        rescale_by_occupancy=_parse_bool(get(opts, "rescale-by-occupancy", string(DEFAULT_RESCALE_BY_OCCUPANCY))),
        downsampled_dir=get(opts, "downsampled-dir", joinpath(Meris.DATADIR, "downsampled")),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_downsampled_omega_afd(; _parse_cli(ARGS)...)
end
