module Utils

using DataFrames, CairoMakie, JLD2
using FHist, Statistics, SparseArrays, LinearAlgebra
using StatsBase, Distributions, SpecialFunctions
using NLsolve

function df_filter!(df::DataFrame; min_samples=1, min_nreads=1)

    # First remove samples with nreads < min_nreads
    filter!(row -> row.nreads >= min_nreads, df)

    # Consider only nevironments with more than min_samples
    grouped = groupby(df, [:env])
    temp = combine(grouped) do sdf
        if length(unique(sdf.sample_id)) >= min_samples
            return sdf
        else
            return DataFrame(
                env    = Float64[],
                species_id  = Float64[],
                sample_id  = Float64[],
                count  = Int64[],
                nreads  = Int64[],
                )
        end
    end

    good_set  = unique(temp.env)
    filter!(row -> row.env in good_set, df)

    return df
end

function get_counts(df; occ = 0.95)
    
    occ = 1 - occ
    
    # Transform raw data into matrix of counts and vector of nreads
    species = unique(df.species_id)
    samples = unique(df.sample_id)
    
    S, T = length(species), length(samples)
    sm_groups = groupby(df, :sample_id)
    
    counts = zeros(T, S)
    nreads = zeros(T)
    otu_index = Dict(sp => i for (i, sp) in enumerate(species))
    run_index = Dict(sm => i for (i, sm) in enumerate(samples))
    
    for g in sm_groups
        sm = g.sample_id[1]
        i = run_index[sm]
        nreads[i] = g.nreads[1]
        for (sp, val) in zip(g.species_id, g.count)
            j = otu_index[sp]
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
    if isnothing(max_idx)
        max_idx = S
    end
    counts = counts[:, 1:max_idx]

    return counts, nreads
end

function get_frequencies(df; occ = 0.9, rescale=true)

    occ = 1 - occ
    dff = copy(df)
    dff.f = dff.count ./ dff.nreads
    
    # Transform raw data into matrix of counts and vector of nreads
    species = unique(df.species_id)
    samples = unique(df.sample_id)
    
    S, T = length(species), length(samples)
    sm_groups = groupby(dff, :sample_id)
    
    freqs = zeros(T, S)
    otu_index = Dict(sp => i for (i, sp) in enumerate(species))
    run_index = Dict(sm => i for (i, sm) in enumerate(samples))
    
    for g in sm_groups
        sm = g.sample_id[1]
        i = run_index[sm]
        for (sp, val) in zip(g.species_id, g.f)
            j = otu_index[sp]
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
    if isnothing(max_idx)
        max_idx = S
    end
    freqs = freqs[:, 1:max_idx]

    if rescale # Multiply by the occupancy
        zero_counts = sum(freqs .!= 0, dims=1)
        freqs  .*= zero_counts ./ T
    end

    return freqs
end

function check_occupancy_thresh(df; occ=0.95)
    sam_ids = unique(df.sample_id)
    T = length(sam_ids)
    ct_map = countmap(df.species_id)
    occ = []
    
    for key in keys(ct_map)
        push!(occ, ct_map[key] / T)
    end

    return length(occ)
end

function reads_distribution(df; bins=30, plot=false, verbose=false, save=false, filename="temp")

    envs = unique(df.env)
    reads_out = Dict()

    for env in envs
        if verbose
            println(env)
        end
        sub   = df[df.env .== env, :]
        
        data = []
        group = groupby(sub, :sample_id)
        for g in group
            push!(data, g.nreads[1])
        end
        
        data = log.(data)
        data = (data .- mean(data)) ./ std(data)
        
        bmin, bmax = extrema(data)
        Δb = (bmax - bmin) / 30
        fh = FHist.Hist1D(data, binedges=bmin:Δb:bmax)
        ctrs = bincenters(fh)
        pdf  = bincounts(fh) ./ (integral(fh) * Δb)
        
        mask = pdf .> 0
        reads_out[env] = (ctrs[mask], pdf[mask])
    end

    if save
        @save "$filename.jld2" reads_out
    end

    if plot
        fig = Figure(figsize=(900,500))
        ax  = Axis(fig[1, 1]; yscale=log10,
                   xlabel = "z", ylabel = "pdf", 
                   title = "Reads Distribution")
    
        for key in keys(reads_out)
            x, y = afd_out[key]
            sc = scatter!(ax, x, 10 .^ log.(y),
                        label=key,
                        markersize=15,
                        strokewidth = 0.8,
                        strokecolor = :black)
        end

        return reads_out, fig
    end

    return reads_out
end

function compute_AFD(df; occ=0.99, bins=30, plot=false, verbose=false, save=false, filename="temp")
    """
    
    Compute the aggregated frequency distribution (AFD) of occurrences in the DataFrame `df`, grouped by the environment column.
    
    # Arguments
    - `df::DataFrame`: Input DataFrame. Must contain a column named `env` and the data columns for frequency calculation.
    - `occ::Float64=0.99`: Occupancy threshold passed to `get_frequencies` to filter rare events.
    - `bins::Integer=30`: Number of histogram bins for estimating the PDF of z-scores.
    - `plot::Bool=false`: If `true`, returns a Makie `Figure` object along with the AFD data.
    - `verbose::Bool=false`: If `true`, prints the name of each environment as it is processed.
    - `save::Bool=false`: If `true`, saves the resulting AFD dictionary to a JLD2 file named `"<filename>.jld2"`.
    - `filename::String="temp"`: Base name for the output file if `save=true`.
    
    # Returns
    - `afd_out::Dict{Any, Tuple{Vector{Float64}, Vector{Float64}}}`: A dictionary mapping each environment to a tuple `(centers, pdf)`:
        - `centers::Vector{Float64}`: Bin centers of the histogram for non-NaN z-scores.
        - `pdf::Vector{Float64}`: Estimated probability density values at each center.
    - If `plot=true`, returns a tuple `(afd_out, fig)` where `fig` is a Makie `Figure` displaying the AFD curves for each environment.
    """
    
    envs = unique(df.env)
    afd_out = Dict()

    for env in envs
        if verbose
            println(env)
        end
        sub   = df[df.env .== env, :]
        freqs = get_frequencies(sub, occ=occ)

        log_non_zero = [log.(col[col .> 0]) for col in eachcol(freqs)]
        μ = mean.(log_non_zero)
        σ = std.(log_non_zero)
        allz = vcat( [(x .- μ[j]) ./ σ[j] for (j, x) in enumerate(log_non_zero)]... )
        allz = allz[.!isnan.(allz)]

        bmin, bmax = round(minimum(allz)), round(maximum(allz))
        Δb = (bmax - bmin) / bins
        fh = FHist.Hist1D(allz, binedges=bmin:Δb:bmax)
        ctrs = bincenters(fh)
        pdf  = bincounts(fh) ./ (integral(fh) * Δb)

        mask = pdf .> 0
        afd_out[env] = (ctrs[mask], pdf[mask])
    end

    if save
        @save "$filename.jld2" afd_out
    end

    if plot
        fig = Figure(figsize=(900,500))
        ax  = Axis(fig[1, 1]; yscale=log10,
                   xlabel = "z", ylabel = "pdf", 
                   title = "AFD")
    
        for key in keys(afd_out)
            x, y = afd_out[key]
            sc = scatter!(ax, x, 10 .^ log.(y),
                        label=key,
                        markersize=15,
                        strokewidth = 0.8,
                        strokecolor = :black)
        end

        leg = Legend(fig, ax; orientation = :vertical)
        fig[1, 2] = leg 

        return afd_out, fig
    end

    return afd_out
end

function compute_TL(df; occ=0.99, bins=30, plot=false, verbose=false, save=false, filename="temp", interpretation=1)
    """
    This function computes Taylor's Law from a given DataFrame `df`, which contains species count data across different environments. The function performs data processing, aggregation, and optionally plots and saves the results.
    
    ### Arguments
    - `df::DataFrame`: The input data frame containing the species count data. The `df` should have a column `env` representing different environments, and other columns representing species counts.
    - `occ::Float64`: The occupancy threshold (default: `0.99`). This parameter is used to filter the data for species that appear in a given proportion of samples.
    - `bins::Int`: The number of bins for binning the mean values (default: `30`). The function divides the x-axis (log-transformed means) into `bins` number of intervals.
    - `plot::Bool`: A flag to indicate whether to generate a plot of Taylor's Law (default: `false`).
    - `verbose::Bool`: A flag for printing progress information (default: `false`).
    - `save::Bool`: A flag to save the output to a file (default: `false`).
    - `filename::String`: The filename to save the output if `save=true` (default: `"temp"`).
    
    ### Returns
    - If `plot=false`: Returns a dictionary `taylor_out`, where the keys are environment labels and the values are tuples of x and y values representing the log-transformed means and variances for each environment.
    - If `plot=true`: Returns both the dictionary `taylor_out` and the figure object `fig` for plotting the results.
    """
    
    envs = unique(df.env)
    taylor_out = Dict()

    for env in envs
        if verbose
            println(env)
        end
        sub   = df[df.env .== env, :]
        counts, nreads = get_counts(sub, occ=occ)
        T, S = size(counts)
        
        if S < 4
            continue
        end

        # Compute mean and var for each species
        if interpretation == 1
            mean_data = sum(counts ./ nreads, dims=1) ./ T
            var_data = sum(counts .* (counts .- 1) ./ (nreads .* (nreads .- 1)), dims=1) ./ T .- (sum(counts ./ nreads, dims=1) ./ T) .^ 2
            mask = var_data .> 0
    
            if length(var_data[mask]) < 2
                continue
            end
        elseif interpretation == 2
            mean_data = [mean(x) for x in eachcol(counts ./ nreads)]
            var_data = [var(x) for x in eachcol(counts ./ nreads)]
            mask = var_data .> 0
        end
    
        # Log transform: it's easier to fit power laws in log-space
        log_mean = log.(mean_data[mask])
        log_var = log.(var_data[mask])
    
        # Bin x-axis (means) and aggregate y-axis (variances)
        bmin = minimum(log_mean)
        bmax = maximum(log_mean)
        Δb = (bmax - bmin) / bins
        binedges = bmin:Δb:bmax
        xx = 0.5 .* (binedges[2:end] .+ binedges[1:end-1])
        yy = [mean(log_var[(log_mean .>= binedges[i]) .& (log_mean .< binedges[i+1])]) for i in 1:length(binedges)-1]

        taylor_out[env] = (xx, yy)
    end

    if save
        @save "$filename.jld2" taylor_out
    end

    if plot
        fig = Figure(figsize=(900,500))
        ax  = Axis(fig[1, 1];
                   xlabel = "log(μ)", ylabel = "log(σ^2)", 
                   title = "Taylor's Law")
    
        for key in keys(taylor_out)
            x, y = taylor_out[key]
            sc = scatter!(ax, x, y,
                        label=key,
                        markersize=15,
                        strokewidth = 0.8,
                        strokecolor = :black)
        end

        leg = Legend(fig, ax; orientation = :vertical)
        fig[1, 2] = leg 

        return taylor_out, fig
    end

    return taylor_out
end

function compute_MAD(df; c=exp(-10), bins=30, plot=false, verbose=false, save=false, filename="temp")
    """
    Compute the median absolute deviation (MAD) distribution of log‐transformed mean frequencies per environment, optionally saving and/or plotting the results.
    
    # Arguments
    - `df::DataFrame`: The input data frame containing the species count data. The `df` should have a column `env` representing different environments, and other columns representing species counts.
    - `c::Union{Float64, Dict}`: Cutoff parameter: if a `Float64`, the same cutoff is applied to all environments; if a `Dict`, a separate cutoff per environment (accessed as `c[env]`).
    - `bins::Int`: The number of bins for binning the mean values (default: `30`). The function divides the x-axis (log-transformed means) into `bins` number of intervals.
    - `plot::Bool`: A flag to indicate whether to generate a plot of Taylor's Law (default: `false`).
    - `verbose::Bool`: A flag for printing progress information (default: `false`).
    - `save::Bool`: A flag to save the output to a file (default: `false`).
    - `filename::String`: The filename to save the output if `save=true` (default: `"temp"`).
    
    # Returns
    - `mad_out::Dict{Any, Tuple{Vector{Float64}, Vector{Float64}}}`: A dictionary mapping each environment to a tuple `(centers, pdf)`:
        - `centers::Vector{Float64}`: Bin centers of the histogram for non-NaN z-scores.
        - `pdf::Vector{Float64}`: Estimated probability density values at each center.
    - If `plot=true`, returns a tuple `(mad_out, fig)` where `fig` is a Makie `Figure` displaying the MAD curves for each environment.
    """

    envs = unique(df.env)
    mad_out = Dict()
    c_dict = Dict()

    for env in envs
        if verbose
            println(env)
        end
        
        if c isa AbstractFloat
            c_dict[env] = c
        elseif c isa AbstractDict
            c_dict[env] = c[env]
        end
        
        sub   = df[df.env .== env, :]
        counts, nreads = get_counts(sub, occ=0)
        T, S = size(counts)
        
        if S < 4
            continue
        end

        mean_data = [mean(x) for x in eachcol(counts ./ nreads)]
        mean_logs = log.(mean_data)
        
        if !isnothing(c)
            mean_logs = mean_logs[mean_logs .> log(c_dict[env])]
        end

        bmin, bmax = round(minimum(mean_logs)), round(maximum(mean_logs))
        Δb = (bmax - bmin) / bins
        fh = FHist.Hist1D(mean_logs, binedges=bmin:Δb:bmax)

        if !isnothing(c)
            m1 = mean(mean_logs)
            m2 = mean(mean_logs .^ 2)
            μ, σ = compute_MAD_params(m1, m2, c_dict[env])
            
            ctrs = bincenters(fh)
            ctrs .-= μ
            ctrs ./= σ
                
            pdf = bincounts(fh) ./ (integral(fh) * Δb)
            valid = pdf .> 0.0
            erfc_arg = (log(c_dict[env]) - μ) / sqrt(2 * σ^2)
            pdf = pdf[valid] .* (erfc(erfc_arg) / 2) .* σ
        else
            μ, σ = mean(fh), std(fh)
            
            ctrs = bincenters(fh)
            ctrs .-= μ
            ctrs ./= σ
                
            pdf = bincounts(fh) ./ (integral(fh) * Δb)
            valid = pdf .> 0.0
            pdf = pdf[valid] .* sqrt(2 * π * σ)
        end

        mad_out[env] = (ctrs[valid], pdf)
    end

    if save
        @save "$filename.jld2" mad_out
    end

    if plot
        fig = Figure(figsize=(900,500))
        ax  = Axis(fig[1, 1]; yscale=log10,
                   xlabel = "rescaled z", ylabel = "pdf", 
                   title = "MAD")
    
        for key in keys(mad_out)
            x, y = mad_out[key]
            sc = scatter!(ax, x, 10 .^ log.(y),
                        label=key,
                        markersize=15,
                        strokewidth = 0.8,
                        strokecolor = :black)
        end

        leg = Legend(fig, ax; orientation = :vertical)
        fig[1, 2] = leg 

        return mad_out, fig
    end

    return mad_out
end

function compute_pearson(df; occ=0.99, bins=30, plot=false, verbose=false, save=false, filename="temp", interpretation = 1)
    """
    Compute the distribution of Pearson correlations in the DataFrame `df`, grouped by the environment column.
    
    # Arguments
    - `df::DataFrame`: Input DataFrame. Must contain a column named `env` and the data columns for frequency calculation.
    - `occ::Float64=0.99`: Occupancy threshold passed to `get_frequencies` to filter rare events.
    - `bins::Integer=30`: Number of histogram bins for estimating the PDF of z-scores.
    - `plot::Bool=false`: If `true`, returns a Makie `Figure` object along with the AFD data.
    - `verbose::Bool=false`: If `true`, prints the name of each environment as it is processed.
    - `save::Bool=false`: If `true`, saves the resulting AFD dictionary to a JLD2 file named `"<filename>.jld2"`.
    - `filename::String="temp"`: Base name for the output file if `save=true`.
    
    # Returns
    - `corr_out::Dict{Any, Tuple{Vector{Float64}, Vector{Float64}}}`: A dictionary mapping each environment to a tuple `(centers, pdf)`:
        - `centers::Vector{Float64}`: Bin centers of the histogram for non-NaN correlations.
        - `pdf::Vector{Float64}`: Estimated probability density values at each center.
    - If `plot=true`, returns a tuple `(corr_out, fig)` where `fig` is a Makie `Figure` displaying the correlation curves for each environment.
    """
    
    envs = unique(df.env)
    corr_out = Dict()

    for env in envs
        if verbose
            println(env)
        end
        sub = df[df.env .== env, :]

        counts, nreads = get_counts(sub; occ=occ)
        S, T = size(counts, 2), size(counts, 1)
        if S < 4
            continue
        end
        
        # Compute mean and var for each species
        if interpretation == 1
            mean_data = sum(counts ./ nreads, dims=1) ./ T
            var_data = sum(counts .* (counts .- 1) ./ (nreads .* (nreads .- 1)), dims=1) ./ T .- mean_data .^ 2
            mask = var_data .> 0
            var_data = var_data[:, mask[1,:]]
            mean_data = mean_data[:, mask[1,:]]
            
            if length(var_data) < 2
                continue
            end
    
            counts = counts[:, mask[1,:]]
    
            allz = [(sum(counts[:,i] .* counts ./ (nreads .* (nreads .- 1)), dims=1) ./ T .- mean_data[i] .* mean_data) ./ sqrt.(var_data[i] .* var_data) for i in 1:size(counts,2)]
            allz = vcat(allz...)
            for i in 1:size(allz, 1)
                allz[i,i] = NaN
            end
            
        elseif interpretation == 2
            allz = pairwise_correlations(counts ./ nreads)
        end
            
        allz = allz[.!isnan.(allz)]
        allz = allz[abs.(allz) .<= 1]
        
        bmin, bmax = minimum(allz), maximum(allz)
        Δb = (bmax - bmin) / bins
        fh = FHist.Hist1D(allz, binedges=bmin:Δb:bmax)
        ctrs = bincenters(fh)
        pdf  = bincounts(fh) ./ (integral(fh) * Δb)

        mask = pdf .> 0
        corr_out[env] = (ctrs[mask], pdf[mask])
    end

    if save
        @save "$filename.jld2" corr_out
    end

    if plot
        fig = Figure(figsize=(900,500))
        ax  = Axis(fig[1, 1]; yscale=log10,
                   xlabel = "ρ", ylabel = "pdf", 
                   title = "Pearson")
    
        for key in keys(corr_out)
            x, y = corr_out[key]
            sc = scatter!(ax, x, 10 .^ log.(y),
                        label=key,
                        markersize=15,
                        strokewidth = 0.8,
                        strokecolor = :black)
        end

        leg = Legend(fig, ax; orientation = :vertical)
        fig[1, 2] = leg 

        return corr_out, fig
    end

    return corr_out
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

"""
symmetric_index(hist; skip=0, normalize=true)

hist = (centers, freqs) or [centers, freqs]
Compares mirrored bins by index (i vs end-i+1).

- skip: ignore that many pairs near the edges
- normalize: divide by total mass so the score is in [0, 2] (≈ scale-free)

Returns: s::Float64
"""
function symmetric(hist; skip::Int=0, normalize::Bool=true)
    centers, freqs = hist
    L = length(freqs)
    @assert L == length(centers) "centers and freqs must match length"
    @assert L ≥ 2 "need at least 2 bins"
    # how many pairs?
    pairs = (L ÷ 2) - skip
    pairs ≤ 0 && return 0.0

    # first half vs reversed second half
    left  = @view freqs[1+skip : skip+pairs]
    right = @view freqs[end-skip : -1 : end-skip-pairs+1]

    s = sum(abs.(left .- right))

    if normalize
        # if freqs are counts, divide by total counts; if densities, this is still a scale
        tot = max(sum(freqs), eps())
        s /= tot
    end
    return s
end


function pairwise_correlations(M::AbstractMatrix)
    # 1) compute the full corr matrix (n×n)
    C = cor(M)  

    # 2) extract just the upper‐triangle, i<j
    #    `triu(mask, k=1)` is a Bool mask with ones for i<j
    mask = triu(trues(size(C)), 1)
    return C[mask]
end

#### ZOO OF DISTRIBUTIONS ####
function lrg(z, α)
    return α*sqrt(trigamma(α)) .* z .+ α*digamma(α) .- exp.(z .* sqrt(trigamma(α)) .+ digamma(α)) .+ 0.5*log(trigamma(α)) .- loggamma(α)
end

function lre(z)
    g = 0.57721566490153286060
    m = -g
    s = π ^ 2 / 6 - g ^ 2
    return log(s) .+ s .* z .+ m .- exp.(s .*z .+ m)
end

function lrl(z, b)
    s = sqrt(trigamma(1) + trigamma(b))
    m = digamma(1) - digamma(b)
    return log(s * b) .+ z .* s .+ m .- (b + 1) .* log.(1 .+ exp.(z .* s .+ m))
end

function lrln(z, σ)
    return -z .^ 2 ./ 2 .- log(sqrt(σ^2 * 2 * π))
end
##############################


end # End module