module TLPrediction

using Meris
using DataFrames
using LsqFit
using Random
using Statistics

"""
    TLPrediction

Utilities for Taylor's-law prediction experiments. The module builds
class-specific count tables, estimates component-level Taylor coefficients
from repeated sample partitions, bins those coefficients on a log scale, and
tests whether the training coefficients predict variance in held-out samples.
"""

# Convert a samples-by-components count matrix into a long DataFrame.
function counts_to_df(ds)
    df = DataFrame(ds, :auto)
    df = stack(
        df,
        names(df),
        variable_name = :component_id,
        value_name = :counts
    )

    df.sample_id = repeat(1:size(ds, 1), outer=size(ds, 2))
    return df
end

# Split a count DataFrame into two sample-disjoint subsets.
function split_df_by_sample(df::DataFrame; rng=Random.default_rng(), frac=0.5)
    samples = shuffle(rng, unique(df.sample_id))
    n1 = floor(Int, frac * length(samples))

    samples1 = Set(samples[1:n1])
    mask1 = in.(df.sample_id, Ref(samples1))

    df1 = df[mask1, :]
    df2 = df[.!mask1, :]

    return df1, df2
end

# Summarize each component by mean, variance, observation count, and Taylor coefficient.
function summary_df(df::DataFrame)
    sdf = combine(
        groupby(df, :component_id),
        :counts => mean => :mean,
        :counts => var => :var,
        :counts => length => :nobs
    )

    transform!(
        sdf,
        [:mean, :var] => ByRow((m, v) -> m > 0 ? (v - m) / m^2 : NaN) => :coeff
    )

    filter!(:coeff => c -> isfinite(c) && c > 1e-5, sdf)

    return sdf
end

# Estimate stable component coefficients across repeated random sample splits.
function repeated_coeff_df(
    df::DataFrame;
    nreps::Int = 100,
    seed::Int = 123,
    frac::Float64 = 0.5,
    min_reps::Int = 10
)
    rng = MersenneTwister(seed)

    all_coeffs = DataFrame()

    for r in 1:nreps
        df1, _ = split_df_by_sample(df; rng=rng, frac=frac)

        sdf = select(summary_df(df1), :component_id, :coeff)
        sdf.rep = fill(r, nrow(sdf))

        append!(all_coeffs, sdf)
    end

    coeff_summary = combine(
        groupby(all_coeffs, :component_id),
        :coeff => mean => :coeff,
        :coeff => std => :coeff_std,
        :coeff => length => :nreps
    )

    filter!(:nreps => n -> n >= min_reps, coeff_summary)

    return coeff_summary, all_coeffs
end

# Bin component coefficients evenly in log10 space and attach bin metadata.
function bin_logcoeff(df::DataFrame, nbins::Int)
    d = copy(df)

    cmin, cmax = extrema(log10.(d.coeff))
    edges = range(cmin, cmax; length=nbins + 1)

    d.coeff_bin = clamp.(searchsortedlast.(Ref(edges), log10.(d.coeff)), 1, nbins)

    d.bin_left_log = edges[d.coeff_bin]
    d.bin_right_log = edges[d.coeff_bin .+ 1]
    d.bin_center_log = (d.bin_left_log .+ d.bin_right_log) ./ 2

    d.bin_left = 10.0 .^ d.bin_left_log
    d.bin_right = 10.0 .^ d.bin_right_log
    d.bin_center = 10.0 .^ d.bin_center_log

    bin_summary = combine(
        groupby(
            d,
            [
                :coeff_bin,
                :bin_left_log,
                :bin_right_log,
                :bin_center_log,
                :bin_left,
                :bin_right,
                :bin_center,
            ]
        ),
        :component_id => length => :ncomponents
    )

    sort!(bin_summary, :coeff_bin)

    component_bins = select(
        d,
        :component_id,
        :coeff,
        :coeff_std,
        :coeff_bin,
    )

    return innerjoin(bin_summary, component_bins, on=:coeff_bin)
end

# Build one count DataFrame per class, optionally downsampling before filtering by occupancy.
function divide_in_class(df::DataFrame; downsample=false, N=nothing, occ=0.05)
    classes = unique(df[:, :class])
    class_dfs = Dict{eltype(classes), DataFrame}()

    for c in classes
        cdf = df[df[:, :class] .== c, :]

        if downsample
            raw_counts, _ = Meris.DataTools.get_counts(cdf; occ=0.0)
            count_totals = vec(sum(raw_counts, dims=2))
            positive_totals = count_totals[count_totals .> 0]
            isempty(positive_totals) && continue

            max_feasible_N = floor(Int, minimum(positive_totals))
            N_c = isnothing(N) ? max_feasible_N : min(Int(N), max_feasible_N)

            _cdf = Meris.DataTools.downsample(cdf; N=N_c, class=c)
            ds_cdf = counts_to_df(Meris.DataTools.order_by_occ(_cdf; occ=occ))
            ds_cdf.nreads .= N_c
        else
            counts, _ = Meris.DataTools.get_counts(cdf; occ=occ)
            ds_cdf = counts_to_df(counts)
        end

        class_dfs[c] = ds_cdf
    end

    return class_dfs
end

# Partition every class-specific DataFrame into train and test sample subsets.
function partition_per_class(class_dfs::Dict; rng=Random.default_rng(), frac=0.5)
    train_dfs = Dict{keytype(class_dfs), DataFrame}()
    test_dfs = Dict{keytype(class_dfs), DataFrame}()

    for (c, cdf) in class_dfs
        train_df, test_df = split_df_by_sample(cdf; rng=rng, frac=frac)
        train_dfs[c] = train_df
        test_dfs[c] = test_df
    end

    return train_dfs, test_dfs
end

# Compute binned training coefficient summaries for each class.
function train_summary(
    class_dfs::Dict;
    nreps::Int = 100,
    seed::Int = 123,
    frac::Float64 = 0.5,
    min_reps::Int = 10,
    nbins::Int = 50
)
    train_summaries = Dict{keytype(class_dfs), DataFrame}()

    for (c, cdf) in class_dfs
        coeff_summary, _ = repeated_coeff_df(
            cdf;
            nreps=nreps,
            seed=seed,
            frac=frac,
            min_reps=min_reps
        )

        train_summaries[c] = bin_logcoeff(coeff_summary, nbins)
    end

    return train_summaries
end

# Join held-out component statistics to training bins and fit each coefficient bin.
function test_summary(test_dfs::Dict, train_summaries::Dict; nbins=100)
    test_summaries = Dict{keytype(test_dfs), DataFrame}()

    for (c, cdf) in test_dfs
        train_summary = train_summaries[c]
        test_summary = select(
            summary_df(cdf),
            [:component_id, :mean, :var]
        )
        component_bins = innerjoin(
            test_summary,
            train_summary,
            on=:component_id
        )

        test_summaries[c] = fit_test(component_bins)
    end

    return test_summaries
end

# Fit var = mean + C * mean^2 within each coefficient bin and compare to C estimates.
function fit_test(test_binned::DataFrame)
    model(x, p) = x .+ p[1] .* x .^ 2

    results = DataFrame(
        coeff_bin = Int[],
        C_est = Float64[],
        C_est_err = Float64[],
        C_fit = Float64[],
        C_fit_err = Float64[]
    )

    for g in groupby(test_binned, :coeff_bin)
        n = nrow(g)

        if n < 10
            continue
        end

        x = g.mean
        y = g.var

        # Least-squares-motivated coefficient average for the bin.
        w = g.mean .^ 4
        C_est, C_est_err = weighted_mean_and_sd(g.coeff, w)

        try
            fit = curve_fit(model, x, y, [C_est])

            C_fit = fit.param[1]
            C_fit_err = stderror(fit)[1]

            push!(
                results,
                (
                    first(g.coeff_bin),
                    C_est,
                    C_est_err,
                    C_fit,
                    C_fit_err
                )
            )
        catch err
            @warn "Fit failed for bin $(first(g.coeff_bin))" err
        end
    end

    return results
end

# Compute a weighted coefficient mean and weighted standard deviation.
function weighted_mean_and_sd(c, w)
    W = sum(w)
    m = sum(w .* c) / W
    s = sqrt(sum(w .* (c .- m) .^ 2) / W)
    return m, s
end

end # End module
