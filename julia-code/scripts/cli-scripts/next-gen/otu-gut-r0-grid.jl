#= Grid analysis for P(R0 > 1) over sigma_B and mean(beta). =#

using DataFrames
using Dates
using Distributions
using JLD2
using LinearAlgebra
using Random
using SparseArrays
using Statistics
using Base.Threads

using Meris

const DEFAULT_CLASS = "GUT1"
const DEFAULT_OUTFILE = joinpath(Meris.DATADIR, "next-gen", "otu-gut1-r0-grid.jld2")
const DEFAULT_BIOMASS_MEAN = "1000.0"
const DEFAULT_BIOMASS_SIGMA_MIN = "0.1"
const DEFAULT_BIOMASS_SIGMA_MAX = "3.0"

function parse_args(args)
    options = Dict(
        "class" => DEFAULT_CLASS,
        "outfile" => DEFAULT_OUTFILE,
        "seed" => "123",
        "connectivity" => "0.01",
        "biomass-mean" => DEFAULT_BIOMASS_MEAN,
        "biomass-sigma-min" => DEFAULT_BIOMASS_SIGMA_MIN,
        "biomass-sigma-max" => DEFAULT_BIOMASS_SIGMA_MAX,
        "biomass-sigmas" => "",
        "biomass-means" => "",
        "beta-mean-min" => "1e-5",
        "beta-mean-max" => "1e-1",
        "beta-means" => "",
        "beta-sigma" => "1.0",
        "diagonal-factor" => "10.0",
        "n-grid" => "20",
        "n-runs" => "10",
        "n-threads" => "10",
        "maxiter" => "1_000",
        "tol" => "1e-8",
        "save-r0" => "false",
    )

    for arg in args
        if arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code --threads=10 julia-code/scripts/cli-scripts/next-gen/otu-gut-r0-grid.jl [options]

            Options:
              --class=NAME              OTU class to analyze. Default: $(DEFAULT_CLASS)
              --outfile=PATH            Output JLD2 file. Default: $(DEFAULT_OUTFILE)
              --seed=N                  Random seed. Default: 123
              --c=P                     Alias for --connectivity.
              --connectivity=P          Off-diagonal Erdos-Renyi edge probability. Default: 0.01
              --biomass-mean=M          Arithmetic mean \\bar{B} of B. Default: $(DEFAULT_BIOMASS_MEAN)
              --biomass-means=a,b,c     Explicit \\bar{B} values, scalar or one per sigma_B.
              --biomass-sigma-min=S     Minimum LogNormal sigma_B. Default: $(DEFAULT_BIOMASS_SIGMA_MIN)
              --biomass-sigma-max=S     Maximum LogNormal sigma_B. Default: $(DEFAULT_BIOMASS_SIGMA_MAX)
              --biomass-sigmas=a,b,c    Explicit sigma_B values; overrides min/max/n-grid.
              --beta-mean-min=M         Minimum arithmetic mean(beta). Default: 1e-5
              --beta-mean-max=M         Maximum arithmetic mean(beta). Default: 1e-1
              --beta-means=a,b,c        Explicit comma-separated mean(beta) values; overrides min/max/n-grid.
              --beta-sigma=S            LogNormal log-space sigma for beta shape. Default: 1
              --diagonal-factor=X       Multiplier applied to beta_ii after drawing beta. Default: 10
              --n-grid=N                Number of grid values per axis. Default: 20
              --n-runs=N                Number of beta matrix runs. Default: 10
              --n-threads=N             Expected Julia thread count. Default: 10
              --maxiter=N               Power iteration limit per sample. Default: 1000
              --tol=X                   Power iteration tolerance. Default: 1e-8
              --save-r0=true|false      Save run/sample/grid-level R0 values. Default: false
            """)
            exit(0)
        elseif startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "="; limit=2)
            key == "c" && (key = "connectivity")
            key == "runs" && (key = "n-runs")
            key == "biomass-sigma" && (key = "biomass-sigmas")
            haskey(options, key) || error("Unknown option: --$key")
            options[key] = value
        else
            error("Unknown argument: $arg")
        end
    end

    parse_optional_grid(value) = isempty(value) ? nothing : parse.(Float64, split(value, ","))

    return (
        class = options["class"],
        outfile = options["outfile"],
        seed = parse(Int, options["seed"]),
        connectivity = parse(Float64, options["connectivity"]),
        biomass_mean = parse(Float64, options["biomass-mean"]),
        biomass_means = parse_optional_grid(options["biomass-means"]),
        biomass_sigma_min = parse(Float64, options["biomass-sigma-min"]),
        biomass_sigma_max = parse(Float64, options["biomass-sigma-max"]),
        biomass_sigmas = parse_optional_grid(options["biomass-sigmas"]),
        beta_mean_min = parse(Float64, options["beta-mean-min"]),
        beta_mean_max = parse(Float64, options["beta-mean-max"]),
        beta_means = parse_optional_grid(options["beta-means"]),
        beta_sigma = parse(Float64, options["beta-sigma"]),
        diagonal_factor = parse(Float64, options["diagonal-factor"]),
        n_grid = parse(Int, replace(options["n-grid"], "_" => "")),
        n_runs = parse(Int, replace(options["n-runs"], "_" => "")),
        n_threads = parse(Int, replace(options["n-threads"], "_" => "")),
        maxiter = parse(Int, replace(options["maxiter"], "_" => "")),
        tol = parse(Float64, options["tol"]),
        save_r0 = lowercase(options["save-r0"]) in ("true", "1", "yes"),
    )
end

function validate_options(options)
    0.0 <= options.connectivity <= 1.0 || error("--connectivity must be between 0 and 1")
    options.biomass_sigma_min >= 0.0 || error("--biomass-sigma-min must be non-negative")
    options.biomass_sigma_max >= options.biomass_sigma_min || error("--biomass-sigma-max must be >= --biomass-sigma-min")
    if options.biomass_sigmas !== nothing
        all(>=(0.0), options.biomass_sigmas) || error("--biomass-sigmas values must be non-negative")
        length(options.biomass_sigmas) > 0 || error("--biomass-sigmas cannot be empty")
    end
    if options.biomass_means !== nothing
        all(>(0.0), options.biomass_means) || error("--biomass-means values must be positive")
        length(options.biomass_means) > 0 || error("--biomass-means cannot be empty")
    end
    options.biomass_mean > 0.0 || error("--biomass-mean must be positive")
    options.beta_mean_min > 0.0 || error("--beta-mean-min must be positive")
    options.beta_mean_max >= options.beta_mean_min || error("--beta-mean-max must be >= --beta-mean-min")
    if options.beta_means !== nothing
        all(>(0.0), options.beta_means) || error("--beta-means values must be positive")
        length(options.beta_means) > 0 || error("--beta-means cannot be empty")
    end
    options.beta_sigma > 0.0 || error("--beta-sigma must be positive")
    options.diagonal_factor >= 0.0 || error("--diagonal-factor must be non-negative")
    options.n_grid > 1 || error("--n-grid must be greater than 1")
    options.n_runs > 0 || error("--n-runs must be positive")
    options.maxiter > 0 || error("--maxiter must be positive")
    options.tol > 0.0 || error("--tol must be positive")

    if Threads.nthreads() != options.n_threads
        @warn "Julia was not launched with the requested thread count" requested=options.n_threads actual=Threads.nthreads()
    end
end

function log_grid(lo::Real, hi::Real, n::Int)
    return exp10.(range(log10(lo), log10(hi), length=n))
end

function linear_grid(lo::Real, hi::Real, n::Int)
    return collect(range(lo, hi, length=n))
end

function lognormal_mu_from_mean_sigma(mean_value::Real, sigma::Real)
    return log(mean_value) - sigma^2 / 2
end

function lognormal_var_from_mean_sigma(mean_value::Real, sigma::Real)
    return (exp(sigma^2) - 1) * mean_value^2
end

function unit_mean_lognormal(sigma::Real)
    return LogNormal(-sigma^2 / 2, sigma)
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

    supports = [findall(>(0.0), @view counts[r, :]) for r in axes(counts, 1)]
    return (; counts, nreads, sample_ids, component_ids, supports)
end

function build_unit_beta_matrix(ncomponents::Int, rng::AbstractRNG; connectivity::Real,
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
                                  v::AbstractVector)
    mul!(tmp, beta_matrix, v)
    @. out = x * tmp
    return out
end

function perron_eigenvalue(beta_matrix::SparseMatrixCSC, x::AbstractVector;
                           tol::Real=1e-8, maxiter::Int=1_000)
    v = fill(1.0 / length(x), length(x))
    y = similar(v)
    tmp = similar(v)
    lambda = 0.0

    for _ in 1:maxiter
        next_generation_product!(y, tmp, beta_matrix, x, v)
        ysum = sum(y)
        ysum > 0.0 || return 0.0

        @. v = y / ysum
        lambda_next = sum(next_generation_product!(y, tmp, beta_matrix, x, v))

        if abs(lambda_next - lambda) <= tol * max(1.0, abs(lambda_next))
            return lambda_next
        end

        lambda = lambda_next
    end

    @warn "Power iteration reached maxiter" maxiter tol
    return lambda
end

function compute_base_rhos(counts::AbstractMatrix, nreads::AbstractVector, supports,
                           beta_matrix::SparseMatrixCSC; tol::Real=1e-8, maxiter::Int=1_000)
    nsamples = size(counts, 1)
    base_rhos = zeros(Float64, nsamples)
    support_sizes = zeros(Int, nsamples)

    Threads.@threads for r in 1:nsamples
        present = supports[r]
        support_sizes[r] = length(present)
        if isempty(present)
            base_rhos[r] = 0.0
            continue
        end

        x = counts[r, present] ./ nreads[r]
        sample_beta = beta_matrix[present, present]
        base_rhos[r] = perron_eigenvalue(sample_beta, x; tol=tol, maxiter=maxiter)
    end

    return base_rhos, support_sizes
end

function update_probability_counts!(gt1_counts, total_counts, biomass_index::Int,
                                    r0_values, beta_means)
    for j in eachindex(beta_means)
        r0_scaled = beta_means[j] .* r0_values
        gt1_counts[biomass_index, j] += count(>(1.0), r0_scaled)
        total_counts[biomass_index, j] += length(r0_values)
    end
end

function result_dataframe(biomass_means, biomass_sigmas, beta_means, gt1_counts, total_counts)
    result = DataFrame(
        biomass_mu = Float64[],
        biomass_sigma = Float64[],
        biomass_mean = Float64[],
        biomass_var = Float64[],
        beta_mean = Float64[],
        n_gt1 = Int[],
        n_total = Int[],
        probability_gt1 = Float64[],
    )

    for i in eachindex(biomass_sigmas), j in eachindex(beta_means)
        total = total_counts[i, j]
        mean_value = biomass_means[i]
        sigma = biomass_sigmas[i]
        push!(
            result,
            (
                lognormal_mu_from_mean_sigma(mean_value, sigma),
                sigma,
                mean_value,
                lognormal_var_from_mean_sigma(mean_value, sigma),
                beta_means[j],
                gt1_counts[i, j],
                total,
                total == 0 ? NaN : gt1_counts[i, j] / total,
            ),
        )
    end

    return result
end

function main(args=ARGS)
    options = parse_args(args)
    validate_options(options)

    rng_master = Random.Xoshiro(options.seed)
    beta_shape_dist = unit_mean_lognormal(options.beta_sigma)
    beta_means = isnothing(options.beta_means) ?
        log_grid(options.beta_mean_min, options.beta_mean_max, options.n_grid) :
        sort(options.beta_means)
    biomass_sigmas = isnothing(options.biomass_sigmas) ?
        linear_grid(options.biomass_sigma_min, options.biomass_sigma_max, options.n_grid) :
        sort(options.biomass_sigmas)
    biomass_means = if options.biomass_means === nothing
        fill(options.biomass_mean, length(biomass_sigmas))
    elseif length(options.biomass_means) == 1
        fill(first(options.biomass_means), length(biomass_sigmas))
    elseif length(options.biomass_means) == length(biomass_sigmas)
        options.biomass_means
    else
        error("--biomass-means must have either one value or one value per sigma_B")
    end
    biomass_mus = [lognormal_mu_from_mean_sigma(mean_value, sigma) for (mean_value, sigma) in zip(biomass_means, biomass_sigmas)]
    biomass_dists = [LogNormal(mu, sigma) for (mu, sigma) in zip(biomass_mus, biomass_sigmas)]

    classdf = load_class_df(options.class)
    data = count_matrix(classdf)
    nsamples = length(data.sample_ids)

    gt1_counts = zeros(Int, length(biomass_sigmas), length(beta_means))
    total_counts = zeros(Int, length(biomass_sigmas), length(beta_means))
    edge_counts = zeros(Int, options.n_runs)
    support_sizes = zeros(Int, nsamples)
    biomass_records = DataFrame(
        run=Int[], sample_id=String[], biomass_mu=Float64[], biomass_sigma=Float64[],
        biomass_mean=Float64[], biomass_var=Float64[], biomass=Float64[]
    )
    r0_records = options.save_r0 ? DataFrame(
        run=Int[], sample_id=String[], biomass_mu=Float64[], biomass_sigma=Float64[],
        biomass_mean=Float64[], biomass_var=Float64[], biomass=Float64[],
        beta_mean=Float64[], R0=Float64[]
    ) : nothing

    for run in 1:options.n_runs
        run_seed = rand(rng_master, UInt)
        beta_rng = Random.Xoshiro(run_seed)
        biomass_rng = Random.Xoshiro(run_seed ⊻ 0x9e3779b97f4a7c15)

        @info "Building unit-mean beta matrix" run n_runs=options.n_runs components=length(data.component_ids) connectivity=options.connectivity
        beta_matrix = build_unit_beta_matrix(
            length(data.component_ids),
            beta_rng;
            connectivity=options.connectivity,
            beta_dist=beta_shape_dist,
            diagonal_factor=options.diagonal_factor,
        )
        edge_counts[run] = nnz(beta_matrix)

        @info "Computing base sample spectral radii" run samples=nsamples edges=edge_counts[run] threads=Threads.nthreads()
        base_rhos, run_support_sizes = compute_base_rhos(
            data.counts,
            data.nreads,
            data.supports,
            beta_matrix;
            tol=options.tol,
            maxiter=options.maxiter,
        )
        support_sizes .= run_support_sizes

        for i in eachindex(biomass_sigmas)
            biomass = rand(biomass_rng, biomass_dists[i], nsamples)
            r0_without_beta_mean = biomass .* base_rhos
            update_probability_counts!(gt1_counts, total_counts, i, r0_without_beta_mean, beta_means)

            append!(
                biomass_records,
                DataFrame(
                    run=fill(run, nsamples),
                    sample_id=data.sample_ids,
                    biomass_mu=fill(biomass_mus[i], nsamples),
                    biomass_sigma=fill(biomass_sigmas[i], nsamples),
                    biomass_mean=fill(biomass_means[i], nsamples),
                    biomass_var=fill(lognormal_var_from_mean_sigma(biomass_means[i], biomass_sigmas[i]), nsamples),
                    biomass=biomass,
                ),
            )

            if options.save_r0
                for j in eachindex(beta_means)
                    append!(
                        r0_records,
                        DataFrame(
                            run=fill(run, nsamples),
                            sample_id=data.sample_ids,
                            biomass_mu=fill(biomass_mus[i], nsamples),
                            biomass_sigma=fill(biomass_sigmas[i], nsamples),
                            biomass_mean=fill(biomass_means[i], nsamples),
                            biomass_var=fill(lognormal_var_from_mean_sigma(biomass_means[i], biomass_sigmas[i]), nsamples),
                            biomass=biomass,
                            beta_mean=fill(beta_means[j], nsamples),
                            R0=beta_means[j] .* r0_without_beta_mean,
                        ),
                    )
                end
            end
        end

        @info "Finished run" run edge_count=edge_counts[run]
        GC.gc()
    end

    result = result_dataframe(biomass_means, biomass_sigmas, beta_means, gt1_counts, total_counts)
    mkpath(dirname(options.outfile))
    created_at = string(now())
    parameters = (
        seed = options.seed,
        class = options.class,
        connectivity = options.connectivity,
        biomass_mean = options.biomass_mean,
        biomass_means = biomass_means,
        biomass_mus = biomass_mus,
        biomass_sigma_min = options.biomass_sigma_min,
        biomass_sigma_max = options.biomass_sigma_max,
        biomass_sigmas = biomass_sigmas,
        beta_mean_min = options.beta_mean_min,
        beta_mean_max = options.beta_mean_max,
        beta_means = beta_means,
        beta_sigma = options.beta_sigma,
        diagonal_factor = options.diagonal_factor,
        n_grid = options.n_grid,
        n_runs = options.n_runs,
        requested_threads = options.n_threads,
        actual_threads = Threads.nthreads(),
        maxiter = options.maxiter,
        tol = options.tol,
        ncomponents = length(data.component_ids),
        nsamples = nsamples,
        edge_counts = edge_counts,
        support_sizes = support_sizes,
    )

    if options.save_r0
        jldsave(
            options.outfile;
            result,
            r0_records,
            biomass_records,
            biomass_mus,
            biomass_means,
            biomass_sigmas,
            beta_means,
            gt1_counts,
            total_counts,
            parameters,
            created_at,
        )
    else
        jldsave(
            options.outfile;
            result,
            biomass_records,
            biomass_mus,
            biomass_means,
            biomass_sigmas,
            beta_means,
            gt1_counts,
            total_counts,
            parameters,
            created_at,
        )
    end

    @info "Saved grid result" file=options.outfile rows=nrow(result)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
