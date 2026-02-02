module Zipf

using Meris
using DataFrames, DataFramesMeta, StatsBase
using JLD2, Random
import Meris.MDistributions as MDist

function compute(
        df;
        relative_counts=false,
        filter=true,
        xmins=nothing,
        save=true,
        filename="zipf.jld2"
    )

    Random.seed!(1234)

    classes = unique(df.class)

    zipf_d = Dict()
    heaps_d = Dict()
    boxplot_d = Dict()
    fit_d = Dict()
    for class in classes
        sdf = df[df.class .== class, :]

        # Compute Heaps'
        heap_df = heaps(sdf)
        samplesize = heap_df.documentsize
        vocabsize = heap_df.vocabularysize
        heaps_d[class] = (samplesize = samplesize, vocabsize = vocabsize)

        # Compute Zipf
        agg_df = aggregate_samples(sdf)
        (relative_counts) && (agg_df.counts ./= agg_df.nreads)
        
        agg_df.ranks .= tiedrank(-agg_df.counts)
        ranks, counts = agg_df.ranks, agg_df.counts
        p = sortperm(ranks)          # permutation that sorts ranks ascending
        ranks  = ranks[p]
        counts = counts[p]

        # Compute parameters for boxplot
        α_vec, β_vec, ε_vec = fit_samples(agg_df; xmins=xmins)
        boxplot_d[class] = α_vec .+ 1

        # Compute TemperedPareto fit
        agg_df.samples_id .= "samp"
        α, β, ε = fit_samples(agg_df; xmins=xmins)
        α, β, ε = α[1], β[1], ε[1]
        fit_d[class] = MDist.TemperedPareto(α, β, ε)
        
        if filter
            mask = counts .> ε
            ranks = ranks[mask]
            counts = counts[mask]
        end

        zipf_d[class] = (ranks = ranks, counts = counts)
    end

    out = (zipf = zipf_d, heaps = heaps_d, boxplot = boxplot_d, fit = fit_d)
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

        fit = MDist.fit(MDist.TemperedPareto, sdf.counts; εs=xmins)
        ε = fit.ε
        push!(ε_vec, ε)
        β = fit.β
        push!(β_vec, β)
        α = fit.α
        push!(α_vec, α)
    end
    
    return (α_vec, β_vec, ε_vec)
end

function heaps(df::DataFrame; sizes=10 .^ collect(1:1e-2:log10(size(df, 1))), rng=Random.Xoshiro(42))
    #/ Construct vocabulary and dictionary
    heaps_df = DataFrame(documentsize=Int[], vocabularysize=Int[])
    component_array = []
    samples = unique(df.sample_id)
    for sample in samples
        sub = df[df.sample_id .== sample, :]
        append!(component_array, sub.component_id)
    end

    shuffle!(rng, component_array)
    for N in Int.(floor.(sizes))
        push!(heaps_df, (N, length(unique(component_array[1:N]))))
    end
    
    return unique(heaps_df)
end

end # End module