module MacroPatterns

using DataFrames, JLD2
using FHist, Statistics
using NLsolve, SpecialFunctions

include("./DataTools.jl")
using .DataTools

"""
    make_hist(data; nbins=20, normalize=true, all_values=false)

Construct a histogram from `data` using `FHist`.

# Arguments
- `data`: Vector of numerical data.
- `nbins`: Number of bins (default = 20).
- `normalize`: If `true`, normalizes the histogram so that the area under the curve is 1.
- `all_values`: If `false`, returns only bins with nonzero counts.

# Returns
A tuple `(centers, pdf)` where:
- `centers`: Vector of bin centers.
- `pdf`: Vector of (possibly normalized) counts or densities.
"""
function make_hist(data; nbins=20, normalize=true, all_values=false)
    bmin, bmax = round(minimum(data)), round(maximum(data))
    Δb = (bmax - bmin) / nbins
    fh = FHist.Hist1D(data, binedges=bmin:Δb:bmax)
    centers = bincenters(fh)
    pdf = bincounts(fh)

    if normalize
        pdf ./= (integral(fh) * Δb)
    end

    if !all_values
        mask = pdf .> 0
        centers = centers[mask]
        pdf = pdf[mask]
    end
    return centers, pdf
end

function compute_AFD(df; occ=0.99, bins=30, verbose=false, save=false, filename="temp")
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
        freqs = DataTools.get_frequencies(sub, occ=occ)

        log_non_zero = [log.(col[col .> 0]) for col in eachcol(freqs)]
        μ = mean.(log_non_zero)
        σ = std.(log_non_zero)
        allz = vcat([(x .- μ[j]) ./ σ[j] for (j, x) in enumerate(log_non_zero)]...)
        allz = allz[.!isnan.(allz)]

        afd_out[env] = make_hist(allz; nbins=bins)
    end

    if save
        @save "$filename.jld2" afd_out
    end

    return afd_out
end

function compute_TL(df; occ=0.99, bins=30, verbose=false, save=false, filename="temp", interpretation=1)
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
        counts, nreads = DataTools.get_counts(sub, occ=occ)
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

    return taylor_out
end

function compute_MAD(df; c=nothing, bins=30, verbose=false, save=false, filename="temp")
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
        counts, nreads = DataTools.get_counts(sub, occ=0)
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

    return mad_out
end

function compute_pearson_distribution(df; occ=0.99, bins=30, verbose=false, save=false, filename="temp", interpretation = 1)
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

        counts, nreads = DataTools.get_counts(sub; occ=occ)
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
        
        corr_out[env] = make_hist(allz; nbins=bins)
    end

    if save
        @save "$filename.jld2" corr_out
    end

    return corr_out
end


### Useful functions
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

function pairwise_correlations(M::AbstractMatrix)
    # 1) compute the full corr matrix (n×n)
    C = cor(M)  

    # 2) extract just the upper‐triangle, i<j
    #    `triu(mask, k=1)` is a Bool mask with ones for i<j
    mask = triu(trues(size(C)), 1)
    return C[mask]
end


end # End module