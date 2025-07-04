module MakeDistributions

using Statistics, StatsBase, Distributions, SpecialFunctions
using NLsolve, LsqFit, FHist
using DataFrames, DataFramesMeta, GLM

function get_frequencies(df; occ = 0.05)

    dff = copy(df)
    dff.f = dff.count ./ dff.nreads
    
    # Transform raw data into matrix of counts and vector of nreads
    otus = unique(df.otu_id)
    runs = unique(df.run_id)
    
    S, T = length(otus), length(runs)
    run_groups = groupby(dff, :run_id)
    
    freqs = zeros(T, S)
    otu_index = Dict(otu => i for (i, otu) in enumerate(otus))
    run_index = Dict(run => i for (i, run) in enumerate(runs))
    
    for g in run_groups
        run = g.run_id[1]
        i = run_index[run]
        for (otu, val) in zip(g.otu_id, g.f)
            j = otu_index[otu]
            freqs[i,j] = val
        end
    end
    
    # Reorder columns counts from most occupied to least occupied
    zero_counts = vec(sum(freqs .== 0, dims=1))
    col_order = sortperm(zero_counts)
    freqs = freqs[:, col_order]
    
    # Filter counts by only consider species with high occupancy
    zero_counts = vec(sum(freqs .== 0, dims=1))
    max_idx = findfirst(>(occ * T), zero_counts)
    freqs = freqs[:, 1:max_idx]

    # Multiply by the occupancy
    zero_counts = sum(freqs .!= 0, dims=1)
    freqs  .*= zero_counts ./ T

    return freqs
end

function get_counts(df; occ = 0.05)
    # Transform raw data into matrix of counts and vector of nreads
    otus = unique(df.otu_id)
    runs = unique(df.run_id)
    
    S, T = length(otus), length(runs)
    run_groups = groupby(df, :run_id)
    
    counts = zeros(T, S)
    nreads = zeros(T)
    otu_index = Dict(otu => i for (i, otu) in enumerate(otus))
    run_index = Dict(run => i for (i, run) in enumerate(runs))
    
    for g in run_groups
        run = g.run_id[1]
        i = run_index[run]
        nreads[i] = g.nreads[1]
        for (otu, val) in zip(g.otu_id, g.count)
            j = otu_index[otu]
            counts[i,j] = val
        end
    end
    
    # Reorder columns counts from most occupied to least occupied
    zero_counts = sum(counts .== 0, dims=1)
    col_order = sortperm(vec(zero_counts))
    counts = counts[:, col_order]
    
    
    # Filter counts by only consider species with high occupancy
    zero_counts = vec(sum(counts .== 0, dims=1))
    max_idx = findfirst(>(occ * T), zero_counts)
    counts = counts[:, 1:max_idx]

    return counts, nreads
end

function first_moment(df; occ=0.05)

    counts, nreads = get_counts(df; occ=occ)
    T = size(counts, 1)
    
    x = copy(counts)
    for i in 1:T
        x[i,:] ./= nreads[i]
    end

    x = sum(x, dims=1) ./ T

    return x
end

function second_moment(df; occ=0.05)

    counts, nreads = get_counts(df; occ=occ)
    T = size(counts, 1)
    
    xx = copy(counts)
    for i in 1:T
        xx[i,:] .*= (xx[i,:] .- 1)
        xx[i,:] ./= nreads[i]*(nreads[i] - 1)
    end

    xx = sum(xx, dims=1) ./ T

    return xx
end

function cross_moments(df; occ=0.05)

    counts, nreads = get_counts(df; occ=occ)
    S, T = size(counts, 2), size(counts, 1)
    
    xx = copy(counts)
    for i in 1:T
        xx[i,:] .*= (xx[i,:] .- 1)
        xx[i,:] ./= nreads[i]*(nreads[i] - 1)
    end
    
    xx = sum(xx, dims=1) ./ T
    
    xy = zeros(size(counts,2),size(counts,2))
    for t in 1:T
        for i in 1:S
            for j in 1:i-1
                xy[i,j] += counts[t,i]*counts[t,j] / (nreads[t] * (nreads[t] - 1))
            end
        end
    end
    
    xy ./= T
    xy .+= xy'

    for i in 1:S
        xy[i,i] = xx[i]
    end

    return xy
end

function get_statistic(df; occ=0.05, cov=false, corr=false)

    x, xx = first_moment(df; occ=occ), second_moment(df; occ=occ)
    
    mean_x = copy(x)
    var_x = xx .- x.^2
    if cov
        xy = cross_moments(df; occ=occ)
        cov_x = xy .- x' * x
        for i in 1:size(cov_x,1)
            cov_x[i,i] = var_x[i]
        end

        return mean_x, var_x, cov_x
    end

    if corr
        corr_x = copy(cov_x)
        for i in 1:size(corr_x,1)
            for j in 1:size(corr_x,1)
                corr_x[i,j] /= sqrt(var_x[i] * var_x[j])
            end
        end

        return mean_x, var_x, cov_x, corr_x
    end

    return mean_x, var_x
end

function make_AFD(data; missing_thresh=size(data, 1), Δb=0.05, env=nothing)

    # Filter data by removing species with low occupancy
    mat = preprocess_matrix(data, make_log=false)
    mask = map(col -> count(ismissing, col) <= missing_thresh, eachcol(mat))
    filtered_mat = mat[:, mask]

    # Remove species that are never present (the previous filtering created empty rows)
    mask = .!map(col -> count(ismissing, col) == size(filtered_mat, 1), eachcol(filtered_mat))
    filtered_mat = filtered_mat[:, mask]
    S = size(filtered_mat, 2)

    # Compute shape parameter β from filtered data:
    # Each species has its own β but according to taylor's law they should all be the same so we take the mean of all βs
    mean_data = [mean(skipmissing(x)) for x in eachcol(filtered_mat)]
    var_data = [var(skipmissing(x)) for x in eachcol(filtered_mat)]
    betas = mean_data.^2 ./ var_data
    mask = .!isnan.(betas) # avois NaNs if any
    β = prod(betas[mask]) ^ (1 / length(betas[mask]))

    # Log-transform and rescaling of filtered data:
    # I want to collapse all data into the same curve so I need to standardize everything
    log_non_zero_data = [log.(skipmissing(x)) for x in eachcol(filtered_mat)]
    log_rescaled_data = [(x .- mean(x)) ./ std(x) for x in log_non_zero_data[mask]]
    log_data = filter(!isnan, vcat(log_rescaled_data...))

    # Compute normalized and recentered Histogram
    bmin = round(minimum(log_data))
    bmax = round(maximum(log_data))
    fh = FHist.Hist1D(log_data, binedges=bmin:Δb:bmax)

    μ = mean(fh)
    σ = std(fh)
    centers = bincenters(fh)
    centers .-= μ
    centers ./= sqrt(2 * σ^2)
    norm_counts = bincounts(fh) ./ (integral(fh) * Δb)

    # Filter non-zero counts
    valid = norm_counts .> 0.0
    yy = log.(norm_counts[valid])
    centers = centers[valid]

    return Dict(
        "hist" => [centers, 10 .^ yy],
        "hparams" => Dict("μ" => μ, "σ" => σ),
        "params" => Dict("β" => β),
        "env" => env
    )
end

function make_Taylor(data; missing_thresh=size(data, 1), Δb=0.05, env=nothing)

    # Filter data by removing species with low occupancy
    mat = preprocess_matrix(data, make_log=false)
    mask = map(col -> count(ismissing, col) <= missing_thresh, eachcol(mat))
    filtered_mat = mat[:, mask]

    # Remove species that are never present (the previous filtering created empty rows)
    mask = .!map(col -> count(ismissing, col) == size(filtered_mat, 1), eachcol(filtered_mat))
    filtered_mat = filtered_mat[:, mask]

    # Compute mean and var for each species
    mean_data = [mean(skipmissing(x)) for x in eachcol(filtered_mat)]
    var_data = [var(skipmissing(x)) for x in eachcol(filtered_mat)]

    # Log transform: it's easier to fit power laws in log-space
    log_mean = log.(mean_data)
    log_var = log.(var_data)

    # Bin x-axis (means) and aggregate y-axis (variances)
    bmin = minimum(log_mean)
    bmax = maximum(log_mean)
    binedges = bmin:Δb:bmax
    centers = 0.5 .* (binedges[2:end] .+ binedges[1:end-1])

    yy = [mean(log_var[(log_mean .>= binedges[i]) .& (log_mean .< binedges[i+1])]) for i in 1:length(binedges)-1]
    centers = centers[isfinite.(yy)]
    yy = yy[isfinite.(yy)]

    # Fit linear model y = αx + q
    func(x, p) = p[1] .* x .+ p[2]
    p0 = [2.0, 0.0]
    fit = curve_fit(func, centers, yy, p0)
    p_fit = fit.param

    return Dict(
        "hist" => [centers, yy],
        "params" => Dict("α" => p_fit[1], "q" => p_fit[2]),
        "env" => env
    )
end

function make_MAD(data; c=exp(-15), missing_thresh=size(data, 1), Δb=0.05, env=nothing)

    # Filter data by removing species with low occupancy
    mat = preprocess_matrix(data, make_log=false)
    mask = map(col -> count(ismissing, col) <= missing_thresh, eachcol(mat))
    filtered_mat = mat[:, mask]

    # Remove species that are never present (the previous filtering created empty rows)
    mask = .!map(col -> count(ismissing, col) == size(filtered_mat, 1), eachcol(filtered_mat))
    filtered_mat = filtered_mat[:, mask]

    # Compute the means and take the log of means bigger than cutoff
    means = [mean(skipmissing(x)) for x in eachcol(filtered_mat)]
    log_data = [log(x) for x in means if x > c]

    # Make the normalized and recenterd Histogram:
    # since we are dealing with a truncated lognormal, we need to compute the moments of the truncated pdf which are different from the moments of the histogram
    bmin = floor(minimum(log_data))
    bmax = ceil(maximum(log_data))
    fh = FHist.Hist1D(log_data, binedges=bmin:Δb:bmax)

    m1 = mean(log_data)
    m2 = mean(log_data .^ 2)
    μ, σ = compute_MAD_params(m1, m2, c)

    centers = bincenters(fh)
    centers .-= μ
    centers ./= sqrt(2 * σ^2)

    norm_counts = bincounts(fh) ./ (integral(fh) * Δb)
    valid = norm_counts .> 0.0
    erfc_arg = (log(c) - μ) / sqrt(2 * σ^2)
    yy = 10 .^ log.((norm_counts[valid]) ./ sqrt(2 / (π * σ^2)) .* erfc(erfc_arg))
    centers = centers[valid]

    return Dict(
        "hist" => [centers, yy],
        "cutoff" => c,
        "hparams" => Dict("μ" => μ, "σ" => σ),
        "env" => env
        )
end

function make_lagCorr(data;
                      missing_thresh = size(data, 1),
                      max_lag = Int(floor(size(data, 1) / 2)),
                      make_log = false, env=nothing)

    mean_corrs, corrs_mat = compute_lagged_autocorrelations(data, max_lag;
                                                             make_log = make_log,
                                                             missing_thresh = missing_thresh)

    return Dict("corrs" => corrs_mat, "mean_corrs" => mean_corrs, "max_lag" => max_lag, "env" => env)
end

function make_lagCrossCorr(data; Δb = 0.01, missing_thresh = size(data, 1), lags = [0], make_log = false, env = nothing)

    corrs = []

    for lag in lags
        cmat = compute_lagged_crosscorrelations(data, lag; make_log = make_log, missing_thresh = missing_thresh)
        push!(corrs, cmat)
    end

    return Dict("cross_corrs" => corrs, "lags" => lags, "env" => env)
end

function make_PSD(data; Δt=1, missing_thresh=size(data, 1), make_log=false, freq_range=nothing, env=nothing)
    mat = preprocess_matrix(data, make_log=make_log)
    if size(mat,1) % 2 != 0
        mat = mat[2:end,:]
    end
    
    N = size(mat, 1)
    N_species = size(data, 2)
    fs = 1 / Δt
    Nf = fs / 2 # Since sample rate is 1 day
    frequencies = (-Int(floor(N/2)):Int(floor(N/2)) - 1) * fs / N # Frequency domain
        
    mean_S = zeros(N) # Initialize array
    otu_count = 0 # Needed for normalization
    for i in 1:N_species
        # Compute non-uniform FFT only for a 'sufficient' number of samples
        if count(ismissing.(mat[:,i])) <= missing_thresh
            otu_count += 1
            
            x = mat[:,i][.!ismissing.(mat[:,i])]
            x .-= mean(x) # Detrend signal to avoid peak at zero frequency and have comparable signals
            
            t_indices = findall(!ismissing, mat[:,i])
            t_normalized = (t_indices .- minimum(t_indices)) ./ N .- 0.5  # The algorithm work for t ∈ [-0.5, 0.5)
            
            p_nfft = NFFT.plan_nfft(t_normalized, N, reltol=1e-9)
            fhat = adjoint(p_nfft) * x
        
            # Compute normalized power spectrum density (periodogram)
            S = abs2.(fhat) .* (Δt / N)
            mean_S .+= S # Mean PSD: this step should give a more accurate result for the PSD supposing that all trajectories are equivalent
        end
    end
        
    # Take only positive frequencies
    positive = frequencies .> 0
    frequencies = frequencies[positive]
    mean_S = mean_S[positive] ./ otu_count

    log_f = log10.(frequencies)
    log_S = log10.(mean_S)

    if !isnothing(freq_range)
        mask = (log_f .>= freq_range[1]) .& (log_f .<= freq_range[2])
        log_f_fit = log_f[mask]
        log_S_fit = log_S[mask]
    else
        log_f_fit = log_f
        log_S_fit = log_S
    end
    
    # Put into a DataFrame and fit linear model: log_S ~ log_k
    plot_df = DataFrame(log_f=log_f_fit, log_S=log_S_fit)
    model = lm(@formula(log_S ~ log_f), plot_df)
    
    # Extract the slope and intercept
    coeffs = coef(model)
    slope = coeffs[2]
    intercept = coeffs[1]

    # standard errors
    se_vec = stderror(model)
    se_intercept = se_vec[1]
    se_slope = se_vec[2]

    return Dict("PSD" => [frequencies, mean_S], "params" => Dict("slope" => (slope, se_slope), "intercept" => (intercept, se_intercept)), "frange" => freq_range, "env" => env)
end

### HELPER FUNCTIONS
# -- Replace 0.0 with `missing` in the data matrix --
function preprocess_matrix(matrix_data::Matrix{Float64}; make_log::Bool=false)
    mat = Matrix{Union{Missing, Float64}}(matrix_data)
    mat[mat .== 0.0] .= missing
    if make_log
        mat = passmissing(log).(mat)
    end
    return mat
end

# -- Custom autocorrelation function that skips missing entries at each lag --
function autocor_skipmissing(x::Union{Vector{Float64}, Vector{Union{Missing, Float64}}}, lag::Int)
    n = length(x)

    if lag == 0
        vals = collect(skipmissing(x))
        return length(vals) > 1 ? cor(vals, vals) : missing
    end

    x1 = x[1:n - lag]
    x2 = x[1 + lag:n]

    valid_pairs = [(x1[i], x2[i]) for i in 1:length(x1) if !ismissing(x1[i]) && !ismissing(x2[i])]

    if length(valid_pairs) < 2
        return missing
    end

    a = first.(valid_pairs)
    b = last.(valid_pairs)

    return cor(a, b)
end

function compute_lagged_autocorrelations(matrix_data::Matrix{Float64}, max_lag::Int64; make_log::Bool=false, missing_thresh::Int64=size(matrix_data, 1))
    mat = preprocess_matrix(matrix_data, make_log=make_log)

    corrs = Vector{Vector{Union{Missing, Float64}}}()

    for i in 1:size(mat, 2)
        x = mat[:, i]
        if count(ismissing.(x)) <= missing_thresh
            c = Vector{Union{Missing, Float64}}()
            for lag in 0:max_lag
                r = autocor_skipmissing(x, lag)
                push!(c, r)
            end
            push!(corrs, c)
        end
    end

    # Sanity check: all vectors same length
    @assert all(length(c) == max_lag + 1 for c in corrs) "Autocorr vectors are not uniform in length"

    # Proper matrix creation
    corrs_mat = hcat(corrs...)

    # Mean across columns (i.e. for each lag)
    mean_corrs = mapslices(x -> mean(skipmissing(x)), corrs_mat; dims=2)

    return mean_corrs, corrs_mat
end

# -- Custom cross correlation function that skips missing entries at each lag --
function crosscor_skipmissing(x::Union{Vector{Float64}, Vector{Union{Missing, Float64}}}, 
                              y::Union{Vector{Float64}, Vector{Union{Missing, Float64}}}, 
                              lag::Int64)
    n = min(length(x), length(y))

    if lag == 0
        vals = [(x[i], y[i]) for i in 1:n if !ismissing(x[i]) && !ismissing(y[i])]
        if length(vals) < 2
            return missing
        end
        a = first.(vals)
        b = last.(vals)
        return cor(a, b)
    elseif lag > 0
        x1 = x[1:n - lag]
        y1 = y[1 + lag:n]
    else
        x1 = x[1 - lag:n]
        y1 = y[1:n + lag]
    end

    valid_pairs = [(x1[i], y1[i]) for i in 1:length(x1) if !ismissing(x1[i]) && !ismissing(y1[i])]

    if length(valid_pairs) < 2
        return missing
    end

    a = first.(valid_pairs)
    b = last.(valid_pairs)

    return cor(a, b)
end

function compute_lagged_crosscorrelations(matrix_data::Matrix{Float64}, lag::Int; make_log::Bool=false, missing_thresh::Int64=size(matrix_data, 1))
    mat = preprocess_matrix(matrix_data, make_log=make_log)
    mask = map(col -> count(ismissing, col) <= missing_thresh, eachcol(mat))
    filtered_mat = mat[:, mask]
    S = size(filtered_mat, 2)
    
    cross_corr = fill(NaN, S, S)  # Use NaN for undefined entries

    for i in 1:S
        for j in 1:i
            if i != j
                x = filtered_mat[:, i]
                y = filtered_mat[:, j]
                val = crosscor_skipmissing(x, y, lag)
                cross_corr[i, j] = ismissing(val) ? NaN : val
                cross_corr[j, i] = ismissing(val) ? NaN : val
            end
        end
    end

    return cross_corr
end

function compute_MAD_params(m1, m2, c)
    function make_system(m1, m2, c)
        return function F!(F, x)
            F[1] = x[1] - m1 + sqrt(2/π) * x[2] * exp(-(log(c) - x[1])^2 / (2 * x[2]^2)) / erfc((log(c) - x[1]) / sqrt(2 * x[2]^2))
            F[2] = x[2]^2 + m1*x[1] + log(c)*m1 - x[1]*log(c) - m2
        end
    end
    
    f! = make_system(m1, m2, c)
    initial_x = [-15.0, 2.0]
    result = nlsolve(f!, initial_x)
    solution = result.zero

    return solution[1], solution[2]
end

end # end module










        
