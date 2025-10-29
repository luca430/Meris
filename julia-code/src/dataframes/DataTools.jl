module DataTools

using DataFrames

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
                env=Float64[],
                species_id=Float64[],
                sample_id=Float64[],
                count=Int64[],
                nreads=Int64[],
            )
        end
    end

    good_set = unique(temp.env)
    filter!(row -> row.env in good_set, df)

    return df
end

function get_counts(df; occ=0.95)

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

function get_frequencies(df; occ=0.9, rescale=true)

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