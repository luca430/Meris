#= Module with tools to parse data =#
module DataTools

using Distributions
using DataFrames, DataFramesMeta
using FHist

#################
### FUNCTIONS ###
"""
Filter standardized DataFrame based on
- minreads: the min. amount of reads within a sample
- minsamples: the min. amount of samples within a class
- mincomponents: the min. amount of distinct components within a class

Note: While in theory we'd want this function to be in-place, it is hard to justify as
      there are some DataFrames functions we with to use that have no in-place equivalent.
"""
function filterdata(
    df::DataFrame;
    minsamples    = 30,
    minreads      = 100_000,
    mincomponents = 100,
    reorder       = true,
    top           = 50
    )
    #~ 1) remove samples with insufficient reads
    df = @subset(df, :nreads .> minreads)
    #~ 2) remove classes with too few distinct components
    class_components = combine(
        groupby(df, :class),
        :component_id => (x -> length(unique(x))) => :ncomponents
    )
    valid_components = class_components.class[class_components.ncomponents .> mincomponents]
    df = semijoin(df, DataFrame(class=valid_components), on=:class)
    #~ 3) remove classes with too few unique samples
    class_samples = combine(
        groupby(df, :class),
        :sample_id => (x -> length(unique(x))) => :nsamples
    )
    valid_samples = class_samples.class[class_samples.nsamples .> minsamples]
    df = semijoin(df, DataFrame(class=valid_samples), on=:class)

    if reorder
        #~ reorder within each class [e.g., a market or a book language]
        sort!(df, [order(:class), order(:nreads, rev=true)])
    end
    #~ select only the top `top` for each :class
    if !isnothing(top)
        summarydf = @chain df begin
            @groupby(:class, :sample_id)
            @combine(:nreads = first(:nreads))
        end
        sort!(summarydf, [order(:class), order(:nreads, rev=true)])
        topdf = combine(groupby(summarydf, :class)) do classdf
            first(classdf, min(top, nrow(classdf)))
        end
        df = semijoin(df, select(topdf, :class, :sample_id), on=[:class, :sample_id])
    end
    #~ Return
    return df
end

"""
Downsample a DataFrame by reducing counts so that each sample has the same `nreads`.
Returns a count matrix with samples in rows and species in columns.
"""
function downsample(df; N=10_000, class=nothing)

    sdf = isnothing(class) ? df : @view df[df.class .== class, :]

    components = unique(sdf.component_id)
    samples = unique(sdf.sample_id)
    S = length(components)

    comp_index = Dict(ci => i for (i, ci) in enumerate(components))
    sample_groups = groupby(sdf, :sample_id)

    nreads = Dict{eltype(samples), Int}()
    occupancy = zeros(Int, S)
    for g in sample_groups
        nreads[g.sample_id[1]] = Int(g.nreads[1])
        for (ci, c) in zip(g.component_id, g.counts)
            if c != 0
                occupancy[comp_index[ci]] += 1
            end
        end
    end

    col_order = sortperm(length(samples) .- occupancy)
    ordered_index = Dict(components[j] => i for (i, j) in enumerate(col_order))
    valid_samples = [sm for sm in samples if nreads[sm] >= N]

    ds_counts = zeros(Int, length(valid_samples), S)
    samp_index = Dict(sm => i for (i, sm) in enumerate(valid_samples))

    for g in sample_groups
        sm = g.sample_id[1]
        !haskey(samp_index, sm) && continue

        i = samp_index[sm]
        remaining = N
        rest = sum(g.counts)

        for (ci, c0) in zip(g.component_id, g.counts)
            remaining == 0 && break

            c = Int(c0)
            if c == 0
                continue
            end

            max_possible = min(c, remaining)
            min_possible = max(0, remaining - (rest - c))

            if min_possible == max_possible
                x = min_possible
            else
                x = rand(Hypergeometric(c, rest - c, remaining))
            end

            ds_counts[i, ordered_index[ci]] = x
            remaining -= x
            rest -= c
        end
    end

    return ds_counts
end

"""
Downsample a standardized DataFrame and return a standardized long DataFrame.
If `N` is not provided, each class is downsampled to its smallest sample `nreads`.
"""
function downsample_df(df::DataFrame; N=nothing, class=nothing)
    if isnothing(class)
        dfs = DataFrame[]
        for c in unique(df.class)
            cdf = downsample_df(df; N=N, class=c)
            nrow(cdf) > 0 && push!(dfs, cdf)
        end
        return isempty(dfs) ? similar(df, 0) : vcat(dfs...)
    end

    sdf = @view df[df.class .== class, :]
    isempty(sdf.class) && return similar(df, 0)

    sample_groups = groupby(sdf, :sample_id)
    sample_nreads = Dict(g.sample_id[1] => Int(g.nreads[1]) for g in sample_groups)
    sample_totals = Dict(g.sample_id[1] => Int(sum(g.counts)) for g in sample_groups)
    feasible_nreads = [
        min(sample_nreads[sm], sample_totals[sm])
        for sm in keys(sample_nreads)
        if min(sample_nreads[sm], sample_totals[sm]) > 0
    ]
    isempty(feasible_nreads) && return similar(df, 0)

    N_c = isnothing(N) ? minimum(feasible_nreads) : min(Int(N), minimum(feasible_nreads))
    N_c <= 0 && return similar(df, 0)

    out = DataFrame(
        class=eltype(sdf.class)[],
        sample_id=eltype(sdf.sample_id)[],
        component_id=eltype(sdf.component_id)[],
        counts=Int[],
        nreads=Int[],
    )

    for g in sample_groups
        min(sample_nreads[g.sample_id[1]], sample_totals[g.sample_id[1]]) < N_c && continue

        remaining = N_c
        rest = sum(g.counts)

        for (component_id, c0) in zip(g.component_id, g.counts)
            remaining == 0 && break

            c = Int(c0)
            if c == 0
                continue
            end

            max_possible = min(c, remaining)
            min_possible = max(0, remaining - (rest - c))

            if min_possible == max_possible
                x = min_possible
            else
                x = rand(Hypergeometric(c, rest - c, remaining))
            end

            if x > 0
                push!(out, (class, g.sample_id[1], component_id, x, N_c))
            end

            remaining -= x
            rest -= c
        end
    end

    return out
end

"""
Order columns of a count matrix by occupancy level.
Optionally, filter by occupancy level.
"""
function order_by_occ(counts; occ=0.0)
    inv_occ = 1 - occ
    # Reorder columns (species) counts from most occupied to least occupied
    zero_counts = sum(counts .== 0, dims=1)
    col_order = sortperm(vec(zero_counts))
    counts = counts[:, col_order]

    # Filter counts by occupancy level
    T, S = size(counts)
    zero_counts = vcat(sum(counts .== 0, dims=1)...)
    max_idx = findfirst(>(inv_occ * T), zero_counts)
    if isnothing(max_idx)
        max_idx = S
    end

    return counts[:, 1:max_idx]
end

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
    bmin, bmax = minimum(data), maximum(data)
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

"Extract a matrix of counts and an array of nreads from standardized DataFrame for a specific level of occupancy."
function get_counts(df; occ=0.95)

    # Transform raw data into matrix of counts and vector of nreads
    components = unique(df.component_id)
    samples = unique(df.sample_id)

    S, T = length(components), length(samples)
    sm_groups = groupby(df, :sample_id)

    counts = zeros(T, S)
    nreads = zeros(T)
    comp_index = Dict(ci => i for (i, ci) in enumerate(components))
    samp_index = Dict(sm => i for (i, sm) in enumerate(samples))

    for g in sm_groups
        sm = g.sample_id[1]
        i = samp_index[sm]
        nreads[i] = g.nreads[1]
        for (ci, val) in zip(g.component_id, g.counts)
            j = comp_index[ci]
            counts[i, j] = val
        end
    end

    counts = order_by_occ(counts; occ=occ)

    return counts, nreads
end

"""
Extract a matrix of frequencies from standardized DataFrame for a specific level of occupancy.
Optionally rescale by occupancy.
"""
function get_frequencies(df; occ=0.95, rescale=false)

    dff = copy(df)
    dff.f = dff.counts ./ dff.nreads

    # Transform raw data into matrix of counts and vector of nreads
    components = unique(df.component_id)
    samples = unique(df.sample_id)

    S, T = length(components), length(samples)
    sm_groups = groupby(dff, :sample_id)

    freqs = zeros(T, S)
    comp_index = Dict(ci => i for (i, ci) in enumerate(components))
    samp_index = Dict(sm => i for (i, sm) in enumerate(samples))

    for g in sm_groups
        sm = g.sample_id[1]
        i = samp_index[sm]
        for (ci, val) in zip(g.component_id, g.f)
            j = comp_index[ci]
            freqs[i, j] = val
        end
    end

    freqs = order_by_occ(freqs; occ=occ)

    if rescale # Multiply by the occupancy
        zero_counts = sum(freqs .!= 0, dims=1)
        freqs .*= zero_counts ./ T
    end

    return freqs
end

end # End module DataTools
