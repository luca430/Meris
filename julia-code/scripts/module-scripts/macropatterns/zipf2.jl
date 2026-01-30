module Zipf

using Meris
using DataFrames, DataFramesMeta, StatsBase
using JLD2, Random
import Meris.MDistributions as MDist

function compute(
        df;
        filter=true,
        xmins=nothing,
        nbins=30,
        save=true,
        filename="zipf.jld2"
    )

    Random.seed!(1234)

    classes = unique(df.class)

    heaps_d = Dict()
    pl_d = Dict()
    cad_d = Dict()
    for class in classes
        sdf = df[df.class .== class, :]
        rdf = deepcopy(sdf)
        rdf.counts .= rdf.counts ./ rdf.nreads
    
        # Compute Heaps' scaling
        hdf = heaps(sdf)
        heaps_d[class] = (N = hdf.documentsize, V = hdf.vocabularysize)
    
        # Compute power law exponents
        α_vec, ε_vec = fit_samples(rdf; xmins=xmins)
    
        # Compute CAD
        agg_df = aggregate_samples(sdf)
        agg_df.counts .= agg_df.counts ./ agg_df.nreads
        agg_df.sample_id .= ""
        α_eff, ε_eff = fit_samples(agg_df; xmins=xmins)
        α_eff, ε_eff = α_eff[1], ε_eff[1]
        pl_d[class] = (α = α_vec, ε = ε_vec, α_eff = α_eff, ε_eff = ε_eff)
        x = filter ? log10.(agg_df.counts[agg_df.counts .> mean(ε_vec)]) : log10.(agg_df.counts)
        h = Meris.DataTools.make_hist(x, nbins=nbins)
        cad_d[class] = (x = h[1], y = h[2])
    end

    out = (heaps = heaps_d, pl = pl_d, cad = cad_d)
    (save) && (@save filename out)
    
    return out
end

#### HELPER ####

function aggregate_samples(df)
    return @chain df begin
        @groupby(:component_id)
        @combine(:sample_id, :agg_counts = sum(:counts), :agg_nreads = sum(:nreads))
        @transform(
            :counts = :agg_counts,
            :nreads = :agg_nreads
        )
    end
end
    
function fit_samples(df; samples_idx=nothing, xmins=nothing)
    
    samples = unique(df.sample_id)
    samples_idx = isnothing(samples_idx) ? collect(1:length(samples)) : samples_idx

    α_vec = []
    β_vec = []
    ε_vec = []
    for sample in samples[samples_idx]
        sdf = df[df.sample_id .== sample, :]

        fit = MDist.fit(MDist.ParetoI, sdf.counts; εs=xmins)
        ε = fit.ε
        push!(ε_vec, ε)
        α = fit.α
        push!(α_vec, α)
    end
    
    return (α_vec, ε_vec)
end

function heaps(
    df::DataFrame;
    sizes = 10 .^ collect(1:1e-2:log10(size(df, 1))),
    nshuffles::Int = 50,
    rng = Random.Xoshiro(42),
)
    # Build component array (typed, not Vector{Any})
    component_array = eltype(df.component_id)[]
    for sample in unique(df.sample_id)
        sub = df[df.sample_id .== sample, :]
        append!(component_array, sub.component_id)
    end

    L = length(component_array)

    # N grid (unique + sorted + valid)
    Ns = sort(unique(Int.(floor.(sizes))))
    Ns = Ns[(Ns .>= 1) .& (Ns .<= L)]
    M = length(Ns)

    sums = zeros(Float64, M)

    # Work buffers
    seen = Set{eltype(component_array)}()

    for _ in 1:nshuffles
        shuffle!(rng, component_array)

        empty!(seen)
        j = 1  # index over Ns
        for i in 1:L
            push!(seen, component_array[i])
            while j <= M && i == Ns[j]
                sums[j] += length(seen)
                j += 1
            end
            j > M && break
        end
    end

    means = sums ./ nshuffles

    return DataFrame(documentsize = Ns, vocabularysize = means)
end

end # End module