module DataTools

using DataFrames

#################
### FUNCTIONS ###
"Filter standardized DataFrame based on min_samples and min_nreads"
function df_filter!(df::DataFrame; min_samples=1, min_nreads=1)

    # First remove samples with nreads < min_nreads
    filter!(row -> row.nreads >= min_nreads, df)

    # Consider only nevironments with more than min_samples
    select!(df, [:class, :component_id, :sample_id, :counts, :nreads])
    grouped = groupby(df, [:class])
    temp = combine(grouped) do sdf
        if length(unique(sdf.sample_id)) >= min_samples
            return sdf
        else
            return DataFrame(
                class=Float64[],
                component_id=Float64[],
                sample_id=Float64[],
                counts=Int64[],
                nreads=Int64[],
            )
        end
    end

    good_set = unique(temp.class)
    filter!(row -> row.class in good_set, df)

    return df
end

"Extract a matrix of counts and an array of nreads from standardized DataFrame for a specific level of occupancy."
function get_counts(df; occ=0.95)

    occ = 1 - occ

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

"""
Extract a matrix of frequencies from standardized DataFrame for a specific level of occupancy.
Optionally rescale by occupancy.
"""
function get_frequencies(df; occ=0.9, rescale=true)

    occ = 1 - occ
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
        freqs .*= zero_counts ./ T
    end

    return freqs
end

end # End module DataTools