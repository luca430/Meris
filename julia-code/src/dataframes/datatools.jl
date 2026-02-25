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
    top           = 10
    )
    #~ filter out samples with insufficient `nreads`
    df = @subset(df, :nreads .> minreads)
    #~ filter out classes with insufficient samples and component diversity
    summarydf = @chain df begin
        @by(
            :class,
            :nsamples = length(:sample_id),
            :ncomponents = length(unique(:component_id))
        )
        @subset(:nsamples .> minsamples, :ncomponents .> mincomponents)
    end
    df = @subset(df, :class .∈ Ref(summarydf.class))

    if reorder
        #~ reorder within each class [e.g., a market or a book language]
        df = combine(groupby(df, :class)) do marketdf
            sort(marketdf, :nreads, rev=true)
        end
    end
    #~ select only the top `top` for each :class
    if !isnothing(top)
        summarydf = @chain df begin
            @by(
                [:class, :sample_id],
                :nreads = first(:nreads)
            )
            #~ Group by `:class` again, rank w.r.t `nreads` and select `top`
            @groupby(:class)
            @transform(:rank = sortperm(:nreads, rev=true))
            @subset(:rank .<= top) #~ select `top` ranks
            @select(Not(:rank))    #~ omit rank as it's not needed [anymore]
        end
        df = semijoin(df, summarydf, on=[:class, :sample_id])
    end
    #~ Return
    return df
end

"""
Filter standardized DataFrame based on min_samples and min_nreads.
"""
function df_filter(
    df::DataFrame;
    minreads::Int=1,
    mincomponents::Int=1,
    minsamplecomponents::Int=1,
    minsamples::Int=30,
)
    #~ filter data
    sdf = @subset(df, :nreads .> minreads)
    sdf = @chain sdf begin
        @groupby(:class)
        @combine(:sample_id, :component_id, :counts, :nreads, :ncomponents = length(unique(:component_id)))
        @subset(:ncomponents .> mincomponents)
        @groupby(:class, :sample_id)
        @combine(:sample_id, :component_id, :counts, :nreads, :ncomponentspersample = length(unique(:component_id)))
        @subset(:ncomponentspersample .> minsamplecomponents)
        @groupby(:class)
        @combine(:sample_id, :component_id, :counts, :nreads, :nsamples = length(unique(:sample_id)))
        @subset(:nsamples .> minsamples)
    end
    return select!(sdf, :class, :sample_id, :component_id, :counts, :nreads)
end

"""
Downsample a DataFrame by reducing counts so that each sample has the same `nreads`.
Returns a count matrix with samples in rows and species in columns.
"""
function downsample(df; N=10_000, class=nothing)

    sdf = deepcopy(df)
    if !isnothing(class)
        sdf = df[df.class.==class, :]
    end

    # Convert df into counts matrix
    counts, nreads = get_counts(sdf, occ=0.0)
    mask = nreads .>= N
    counts = Int.(counts[mask, :])
    nreads = Int.(nreads[mask])

    nrows, ncols = size(counts)
    ds_counts = zeros(Int, nrows, ncols)

    for q in 1:nrows
        row = counts[q, :]
        total = sum(row)
        remaining = N
        rest = total

        for k in 1:ncols
            if remaining == 0
                break
            end

            c = row[k]
            if c == 0
                continue
            end

            # enforce feasibility
            max_possible = min(c, remaining)
            min_possible = max(0, remaining - (rest - c))

            if min_possible == max_possible
                x = min_possible
            else
                x = rand(Hypergeometric(c, rest - c, remaining))
            end

            ds_counts[q, k] = x
            remaining -= x
            rest -= c
        end
    end

    return ds_counts
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
