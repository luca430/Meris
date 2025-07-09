module Utils

using DataFrames, CairoMakie, JLD2
using FHist, Statistics, SparseArrays

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

function get_frequencies(df; occ = 0.9)

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
    max_idx = findfirst(>=(occ * T), zero_counts)
    if isnothing(max_idx)
        max_idx = T
    end
    freqs = freqs[:, 1:max_idx]

    # Multiply by the occupancy
    zero_counts = sum(freqs .!= 0, dims=1)
    freqs  .*= zero_counts ./ T

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
                   xlabel = "log(z)", ylabel = "pdf", 
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


end