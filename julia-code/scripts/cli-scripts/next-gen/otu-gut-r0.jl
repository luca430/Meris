#= Build next-generation matrices for a selected OTU gut class and save its R0 values. =#

using DataFrames
using Dates
using Distributions
using JLD2
using LinearAlgebra
using Random
using SparseArrays

using Meris

const DEFAULT_CLASS = "GUT1"
const DEFAULT_OUTFILE = joinpath(Meris.DATADIR, "next-gen", "otu-gut1-r0.jld2")

function parse_args(args)
    options = Dict(
        "class" => DEFAULT_CLASS,
        "outfile" => DEFAULT_OUTFILE,
        "seed" => "123",
        "connectivity" => "0.01",
        "beta-mean" => "1e-2",
        "beta-sigma" => "1.0",
        "diagonal-factor" => "10.0",
        "biomass-mean" => string(exp(1.5)),
        "biomass-sigma" => "1.0",
        "n-runs" => "1",
        "maxiter" => "1_000",
        "tol" => "1e-8",
    )

    for arg in args
        if arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code julia-code/scripts/cli-scripts/next-gen/otu-gut-r0.jl [options]

            Options:
              --class=NAME           OTU class to analyze. Default: $(DEFAULT_CLASS)
              --outfile=PATH         Output JLD2 file. Default: $(DEFAULT_OUTFILE)
              --seed=N               Random seed. Default: 123
              --c=P                  Alias for --connectivity.
              --connectivity=P       Off-diagonal Erdos-Renyi edge probability. Default: 0.01
              --beta-mean=M          Arithmetic mean of beta_ij. Default: 1e-2
              --beta-sigma=S         LogNormal log-space sigma for beta_ij. Default: 1
              --diagonal-factor=X    Multiplier applied to beta_ii after drawing beta. Default: 10
              --biomass-mean=M       Arithmetic mean \\bar{B} of B^r. Default: exp(1.5)
              --biomass-sigma=S      LogNormal sigma_B for B^r. Default: 1
              --n-runs=N             Number of runs; each run samples a new beta matrix. Default: 1
              --maxiter=N            Power iteration limit per sample. Default: 1000
              --tol=X                Power iteration tolerance. Default: 1e-8

            Parameters can be comma-separated lists, read in order as distinct distributions.
            Scalar values are reused for all distributions. For example:
              --beta-mean=1e-3,1e-2 --biomass-mean=100,1000 --biomass-sigma=0.5,1 --n-runs=5
            """)
            exit(0)
        elseif startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "="; limit=2)
            key == "c" && (key = "connectivity")
            key == "runs" && (key = "n-runs")
            haskey(options, key) || error("Unknown option: --$key")
            options[key] = value
        else
            error("Unknown argument: $arg")
        end
    end

    parse_float_list(value) = parse.(Float64, split(value, ","))
    parse_int_list(value) = parse.(Int, replace.(split(value, ","), "_" => ""))

    return (
        class = options["class"],
        outfile = options["outfile"],
        seed = parse(Int, options["seed"]),
        connectivity = parse_float_list(options["connectivity"]),
        beta_mean = parse_float_list(options["beta-mean"]),
        beta_sigma = parse_float_list(options["beta-sigma"]),
        diagonal_factor = parse_float_list(options["diagonal-factor"]),
        biomass_mean = parse_float_list(options["biomass-mean"]),
        biomass_sigma = parse_float_list(options["biomass-sigma"]),
        n_runs = parse_int_list(options["n-runs"]),
        maxiter = parse(Int, replace(options["maxiter"], "_" => "")),
        tol = parse(Float64, options["tol"]),
    )
end

function _broadcast_parameter(values, n::Int, name::AbstractString)
    length(values) == n && return values
    length(values) == 1 && return fill(first(values), n)
    error("--$name must have either one value or $n values")
end

function distribution_parameters(options)
    lengths = length.([
        options.connectivity,
        options.beta_mean,
        options.beta_sigma,
        options.diagonal_factor,
        options.biomass_mean,
        options.biomass_sigma,
        options.n_runs,
    ])
    n_distributions = maximum(lengths)

    connectivity = _broadcast_parameter(options.connectivity, n_distributions, "connectivity")
    beta_mean = _broadcast_parameter(options.beta_mean, n_distributions, "beta-mean")
    beta_sigma = _broadcast_parameter(options.beta_sigma, n_distributions, "beta-sigma")
    diagonal_factor = _broadcast_parameter(options.diagonal_factor, n_distributions, "diagonal-factor")
    biomass_mean = _broadcast_parameter(options.biomass_mean, n_distributions, "biomass-mean")
    biomass_sigma = _broadcast_parameter(options.biomass_sigma, n_distributions, "biomass-sigma")
    n_runs = _broadcast_parameter(options.n_runs, n_distributions, "n-runs")

    params_df = DataFrame(
        distribution_id = collect(1:n_distributions),
        connectivity = connectivity,
        beta_mean = beta_mean,
        beta_sigma = beta_sigma,
        beta_mu = fill(NaN, n_distributions),
        diagonal_factor = diagonal_factor,
        biomass_mean = biomass_mean,
        biomass_var = [lognormal_var_from_mean_sigma(mean_value, sigma) for (mean_value, sigma) in zip(biomass_mean, biomass_sigma)],
        biomass_mu = [lognormal_mu_from_mean_sigma(mean_value, sigma) for (mean_value, sigma) in zip(biomass_mean, biomass_sigma)],
        biomass_sigma = biomass_sigma,
        n_runs = n_runs,
    )

    return params_df
end

function validate_options(options, params_df::DataFrame)
    all(x -> 0.0 <= x <= 1.0, params_df.connectivity) || error("--connectivity must be between 0 and 1")
    all(>(0.0), params_df.beta_mean) || error("--beta-mean must be positive")
    all(>(0.0), params_df.beta_sigma) || error("--beta-sigma must be positive")
    all(>=(0.0), params_df.diagonal_factor) || error("--diagonal-factor must be non-negative")
    all(>(0.0), params_df.biomass_mean) || error("--biomass-mean must be positive")
    all(>=(0.0), params_df.biomass_sigma) || error("--biomass-sigma must be non-negative")
    all(>(0), params_df.n_runs) || error("--n-runs must be positive")
    options.maxiter > 0 || error("--maxiter must be positive")
    options.tol > 0.0 || error("--tol must be positive")
end

function beta_distribution(mean_value::Real, sigma::Real)
    mu = log(mean_value) - sigma^2 / 2
    return LogNormal(mu, sigma)
end

function lognormal_mu_from_mean_sigma(mean_value::Real, sigma::Real)
    return log(mean_value) - sigma^2 / 2
end

function lognormal_var_from_mean_sigma(mean_value::Real, sigma::Real)
    return (exp(sigma^2) - 1) * mean_value^2
end

function load_class_df(class_name::AbstractString)
    @info "Loading OTU dataset with default filtering"
    df = Meris.OTULoader.load()
    classdf = df[df.class .== class_name, :]
    nrow(classdf) > 0 || error("No rows found for OTU class $class_name")
    return classdf
end

function count_matrix(df::DataFrame)
    component_ids = sort(unique(df.component_id))
    sample_ids = sort(unique(df.sample_id))
    component_index = Dict(component_id => i for (i, component_id) in enumerate(component_ids))
    sample_index = Dict(sample_id => i for (i, sample_id) in enumerate(sample_ids))

    counts = zeros(Float64, length(sample_ids), length(component_ids))
    nreads = zeros(Float64, length(sample_ids))

    for sampledf in groupby(df, :sample_id)
        sample_id = first(sampledf.sample_id)
        row = sample_index[sample_id]
        nreads[row] = first(sampledf.nreads)

        for (component_id, count) in zip(sampledf.component_id, sampledf.counts)
            counts[row, component_index[component_id]] = count
        end
    end

    return (; counts, nreads, sample_ids, component_ids)
end

function build_beta_matrix(ncomponents::Int, rng::AbstractRNG; connectivity::Real,
                           beta_dist::LogNormal, diagonal_factor::Real)
    adjacency = sprand(rng, ncomponents, ncomponents, connectivity) .> 0.0
    rows, cols, _ = findnz(adjacency)
    offdiag = rows .!= cols

    rows = vcat(rows[offdiag], collect(1:ncomponents))
    cols = vcat(cols[offdiag], collect(1:ncomponents))
    beta_values = vcat(
        rand(rng, beta_dist, count(offdiag)),
        diagonal_factor .* rand(rng, beta_dist, ncomponents),
    )

    return sparse(rows, cols, beta_values, ncomponents, ncomponents)
end

function next_generation_product!(out::AbstractVector, tmp::AbstractVector,
                                  beta_matrix::SparseMatrixCSC, x::AbstractVector,
                                  biomass::Real, v::AbstractVector)
    mul!(tmp, beta_matrix, v)
    @. out = biomass * x * tmp
    return out
end

function perron_eigenvalue(beta_matrix::SparseMatrixCSC, x::AbstractVector, biomass::Real;
                           tol::Real=1e-8, maxiter::Int=1_000)
    v = fill(1.0 / length(x), length(x))
    y = similar(v)
    tmp = similar(v)
    lambda = 0.0

    for _ in 1:maxiter
        next_generation_product!(y, tmp, beta_matrix, x, biomass, v)
        ysum = sum(y)
        ysum > 0.0 || return 0.0

        @. v = y / ysum
        lambda_next = sum(next_generation_product!(y, tmp, beta_matrix, x, biomass, v))

        if abs(lambda_next - lambda) <= tol * max(1.0, abs(lambda_next))
            return lambda_next
        end

        lambda = lambda_next
    end

    @warn "Power iteration reached maxiter" maxiter tol
    return lambda
end

function compute_r0_values(counts::AbstractMatrix, nreads::AbstractVector, beta_matrix,
                           rng::AbstractRNG; biomass_dist::LogNormal,
                           tol::Real=1e-8, maxiter::Int=1_000, run::Int=1)
    nsamples = size(counts, 1)
    biomass = rand(rng, biomass_dist, nsamples)
    r0 = zeros(Float64, nsamples)
    support_sizes = zeros(Int, nsamples)

    for r in 1:nsamples
        present = findall(>(0.0), @view counts[r, :])
        support_sizes[r] = length(present)
        if isempty(present)
            r0[r] = 0.0
            continue
        end

        x = counts[r, present] ./ nreads[r]
        sample_beta = beta_matrix[present, present]
        r0[r] = perron_eigenvalue(sample_beta, x, biomass[r]; tol=tol, maxiter=maxiter)

        if r == 1 || r % 10 == 0 || r == nsamples
            @info "Computed sample R0" run sample=r nsamples support=support_sizes[r] edges=nnz(sample_beta)
        end
    end

    return (; r0, biomass, support_sizes)
end

function main(args=ARGS)
    options = parse_args(args)
    distribution_params = distribution_parameters(options)
    validate_options(options, distribution_params)

    rng = Random.Xoshiro(options.seed)

    classdf = load_class_df(options.class)
    data = count_matrix(classdf)
    r0df = DataFrame(
        distribution_id = Int[],
        run = Int[],
        class = String[],
        sample_id = String[],
        biomass = Float64[],
        nspecies = Int[],
        R0 = Float64[],
    )

    edge_count_dfs = DataFrame[]

    for dist_row in eachrow(distribution_params)
        beta_dist = beta_distribution(dist_row.beta_mean, dist_row.beta_sigma)
        biomass_dist = LogNormal(dist_row.biomass_mu, dist_row.biomass_sigma)
        distribution_params[dist_row.distribution_id, :beta_mu] = params(beta_dist)[1]

        edge_counts = zeros(Int, dist_row.n_runs)
        @info "Starting distribution" distribution_id=dist_row.distribution_id n_runs=dist_row.n_runs connectivity=dist_row.connectivity beta_mean=dist_row.beta_mean biomass_mean=dist_row.biomass_mean biomass_sigma=dist_row.biomass_sigma

        for run in 1:dist_row.n_runs
            @info "Building beta matrix" distribution_id=dist_row.distribution_id run n_runs=dist_row.n_runs components=length(data.component_ids) connectivity=dist_row.connectivity
            beta_matrix = build_beta_matrix(
                length(data.component_ids),
                rng;
                connectivity=dist_row.connectivity,
                beta_dist=beta_dist,
                diagonal_factor=dist_row.diagonal_factor,
            )
            edge_counts[run] = nnz(beta_matrix)

            @info "Computing R0 values" distribution_id=dist_row.distribution_id run samples=length(data.sample_ids) edges=edge_counts[run]
            values = compute_r0_values(
                data.counts,
                data.nreads,
                beta_matrix,
                rng;
                biomass_dist=biomass_dist,
                tol=options.tol,
                maxiter=options.maxiter,
                run=run,
            )

            append!(
                r0df,
                DataFrame(
                    distribution_id = fill(dist_row.distribution_id, length(data.sample_ids)),
                    run = fill(run, length(data.sample_ids)),
                    class = fill(options.class, length(data.sample_ids)),
                    sample_id = data.sample_ids,
                    biomass = values.biomass,
                    nspecies = values.support_sizes,
                    R0 = values.r0,
                ),
            )
            GC.gc()
        end

        push!(
            edge_count_dfs,
            DataFrame(
                distribution_id = fill(dist_row.distribution_id, dist_row.n_runs),
                run = collect(1:dist_row.n_runs),
                edge_count = edge_counts,
            ),
        )
    end

    edge_counts = isempty(edge_count_dfs) ?
        DataFrame(distribution_id=Int[], run=Int[], edge_count=Int[]) :
        vcat(edge_count_dfs...)

    mkpath(dirname(options.outfile))
    created_at = string(now())
    parameters = (
        seed = options.seed,
        n_distributions = nrow(distribution_params),
        maxiter = options.maxiter,
        tol = options.tol,
        ncomponents = length(data.component_ids),
    )

    jldsave(
        options.outfile;
        r0df,
        component_ids = data.component_ids,
        distribution_parameters = distribution_params,
        edge_counts,
        parameters,
        created_at,
    )

    @info "Saved R0 values" file=options.outfile samples=nrow(r0df)
    return r0df
end

main()
