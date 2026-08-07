#= Fit cached candidate curves to next-generation R0 distributions. =#

using DataFrames
using Dates
using Distributions
using JLD2
using Printf
using Random
using Statistics

using Meris

const DEFAULT_INPUT = joinpath(Meris.DATADIR, "next-gen", "otu-gut1-r0.jld2")
const DEFAULT_OUTFILE = joinpath(Meris.DATADIR, "next-gen", "otu-gut1-r0-distribution-fit.jld2")

function parse_args(args)
    options = Dict(
        "input" => DEFAULT_INPUT,
        "outfile" => DEFAULT_OUTFILE,
        "distribution-ids" => "",
        "n-bootstrap" => "250",
        "seed" => "123",
        "n-curve" => "400",
    )

    for arg in args
        if arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code julia-code/scripts/cli-scripts/next-gen/fit-otu-gut-r0-lognormal.jl [options]

            Options:
              --input=PATH              Input R0 JLD2 file. Default: $(DEFAULT_INPUT)
              --outfile=PATH            Output fit JLD2 file. Default: $(DEFAULT_OUTFILE)
              --distribution-ids=a,b     Optional comma-separated subset of distribution IDs.
              --n-bootstrap=N           Parametric bootstrap samples for KS GOF. Default: 250
              --seed=N                  Random seed for bootstrap. Default: 123
              --n-curve=N               Number of saved curve points per distribution. Default: 400
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
        outfile = options["outfile"],
        distribution_ids = distribution_ids,
        n_bootstrap = parse(Int, replace(options["n-bootstrap"], "_" => "")),
        seed = parse(Int, options["seed"]),
        n_curve = parse(Int, replace(options["n-curve"], "_" => "")),
    )
end

function load_r0(path::AbstractString)
    isfile(path) || error("Input file not found: $path")
    data = JLD2.load(path)
    haskey(data, "r0df") || error("Expected key `r0df` in $path")
    return data["r0df"], get(data, "parameters", nothing), get(data, "distribution_parameters", nothing)
end

function prepare_r0df!(r0df::DataFrame)
    if :distribution_id ∉ propertynames(r0df)
        r0df.distribution_id = fill(1, nrow(r0df))
    end
    return r0df
end

function filter_distributions(r0df::DataFrame, distribution_ids)
    distribution_ids === nothing && return r0df
    filtered = r0df[in.(r0df.distribution_id, Ref(distribution_ids)), :]
    nrow(filtered) > 0 || error("No rows found for distribution IDs $(join(distribution_ids, ","))")
    return filtered
end

function fit_lognormal(values::AbstractVector{<:Real})
    logs = log.(values)
    mu = mean(logs)
    sigma = sqrt(mean((logs .- mu) .^ 2))
    sigma > 0.0 || error("Cannot fit a lognormal distribution with zero log-space variance")
    return LogNormal(mu, sigma), mu, sigma
end

function fit_gamma(values::AbstractVector{<:Real})
    dist = fit_mle(Gamma, values)
    alpha, theta = params(dist)
    return dist, alpha, theta
end

function ks_distance(values::AbstractVector{<:Real}, dist::Distribution)
    x = sort(collect(values))
    n = length(x)
    F = cdf.(dist, x)
    dplus = maximum((1:n) ./ n .- F)
    dminus = maximum(F .- (0:n-1) ./ n)
    return max(dplus, dminus)
end

function refit_distribution(model::Symbol, values::AbstractVector{<:Real})
    model == :lognormal && return first(fit_lognormal(values))
    model == :gamma && return first(fit_gamma(values))
    error("Unknown model: $model")
end

function bootstrap_ks_pvalue(model::Symbol, values::AbstractVector{<:Real}, fitted::Distribution, observed_ks::Real,
                             rng::AbstractRNG; n_bootstrap::Int)
    n_bootstrap <= 0 && return NaN

    n = length(values)
    exceedances = 0
    for _ in 1:n_bootstrap
        sample = rand(rng, fitted, n)
        refitted = refit_distribution(model, sample)
        exceedances += ks_distance(sample, refitted) >= observed_ks
    end

    return (exceedances + 1) / (n_bootstrap + 1)
end

function log10_density(dist::Distribution, x::Real)
    r0 = 10.0^x
    return log(10) * r0 * pdf(dist, r0)
end

function fit_curve(distribution_id::Int, model::Symbol, dist::Distribution,
                   log10_min::Real, log10_max::Real, n_curve::Int)
    xs = collect(range(log10_min, log10_max, length=n_curve))

    return DataFrame(
        distribution_id = fill(distribution_id, n_curve),
        model = fill(String(model), n_curve),
        logR0 = xs,
        density = log10_density.(Ref(dist), xs),
    )
end

function candidate_row(distribution_id::Int, model::Symbol, values::AbstractVector{<:Real},
                       dist::Distribution, param1::Real, param2::Real, rng::AbstractRNG;
                       n_bootstrap::Int)
    loglik = sum(logpdf.(dist, values))
    n = length(values)
    nparams = 2
    ks = ks_distance(values, dist)
    pvalue = bootstrap_ks_pvalue(model, values, dist, ks, rng; n_bootstrap=n_bootstrap)

    return (
        distribution_id,
        String(model),
        n,
        param1,
        param2,
        loglik,
        2 * nparams - 2 * loglik,
        log(n) * nparams - 2 * loglik,
        ks,
        pvalue,
        false,
    )
end

function fit_all(r0df::DataFrame, rng::AbstractRNG; n_bootstrap::Int=250, n_curve::Int=400)
    fitdf = DataFrame(
        distribution_id = Int[],
        model = String[],
        n = Int[],
        param1 = Float64[],
        param2 = Float64[],
        loglik = Float64[],
        aic = Float64[],
        bic = Float64[],
        ks_distance = Float64[],
        ks_bootstrap_pvalue = Float64[],
        best = Bool[],
    )
    curve_dfs = DataFrame[]

    global_log10 = log10.(r0df.R0)
    log10_min, log10_max = extrema(global_log10)

    for distribution_id in sort(unique(r0df.distribution_id))
        sdf = r0df[r0df.distribution_id .== distribution_id, :]
        values = collect(sdf.R0)

        lognormal_dist, lognormal_mu, lognormal_sigma = fit_lognormal(values)
        gamma_dist, gamma_alpha, gamma_theta = fit_gamma(values)

        lognormal_row = candidate_row(
            distribution_id,
            :lognormal,
            values,
            lognormal_dist,
            lognormal_mu,
            lognormal_sigma,
            rng;
            n_bootstrap=n_bootstrap,
        )
        gamma_row = candidate_row(
            distribution_id,
            :gamma,
            values,
            gamma_dist,
            gamma_alpha,
            gamma_theta,
            rng;
            n_bootstrap=n_bootstrap,
        )

        first_row = nrow(fitdf) + 1
        push!(fitdf, lognormal_row)
        push!(fitdf, gamma_row)

        best_is_lognormal = lognormal_row[7] <= gamma_row[7]
        best_row = best_is_lognormal ? first_row : first_row + 1
        fitdf[best_row, :best] = true

        best_model = best_is_lognormal ? :lognormal : :gamma
        best_dist = best_is_lognormal ? lognormal_dist : gamma_dist
        push!(curve_dfs, fit_curve(distribution_id, best_model, best_dist, log10_min, log10_max, n_curve))
    end

    curve_df = isempty(curve_dfs) ?
        DataFrame(distribution_id=Int[], model=String[], logR0=Float64[], density=Float64[]) :
        vcat(curve_dfs...)

    return fitdf, curve_df
end

function print_gof(fitdf::DataFrame)
    println("Distribution goodness of fit by R0 distribution")
    println("distribution_id,model,n,param1,param2,KS_D,KS_bootstrap_pvalue,AIC,BIC,best")
    for row in eachrow(fitdf)
        @printf(
            "%d,%s,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%s\n",
            row.distribution_id,
            row.model,
            row.n,
            row.param1,
            row.param2,
            row.ks_distance,
            row.ks_bootstrap_pvalue,
            row.aic,
            row.bic,
            row.best,
        )
    end
end

function main(args=ARGS)
    options = parse_args(args)
    options.n_bootstrap >= 0 || error("--n-bootstrap must be non-negative")
    options.n_curve >= 2 || error("--n-curve must be at least 2")

    r0df, parameters, distribution_parameters = load_r0(options.input)
    prepare_r0df!(r0df)
    r0df = filter_distributions(r0df, options.distribution_ids)
    r0df = r0df[r0df.R0 .> 0.0, :]
    nrow(r0df) > 0 || error("No positive R0 values found in $(options.input)")

    rng = Random.Xoshiro(options.seed)
    fitdf, curve_df = fit_all(r0df, rng; n_bootstrap=options.n_bootstrap, n_curve=options.n_curve)

    mkpath(dirname(options.outfile))
    created_at = string(now())
    fit_parameters = (
        input = options.input,
        seed = options.seed,
        n_bootstrap = options.n_bootstrap,
        n_curve = options.n_curve,
    )

    jldsave(
        options.outfile;
        fitdf,
        curve_df,
        parameters,
        distribution_parameters,
        fit_parameters,
        created_at,
    )

    @info "Saved R0 distribution fits" file=options.outfile candidates=nrow(fitdf)
    print_gof(fitdf)
    return fitdf
end

main()
