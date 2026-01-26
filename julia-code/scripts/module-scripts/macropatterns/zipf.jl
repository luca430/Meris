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

    zipf = Dict()
    heaps = Dict()
    boxplot = Dict()
    fit = Dict()
    for class in classes
        sdf = df[df.class .== class, :]

        # Compute Heaps'
        heap_df = computevocabsize(sdf)
        samplesize = heap_df.documentsize
        vocabsize = heap_df.vocabularysize
        heaps[class] = (samplesize = samplesize, vocabsize = vocabsize)

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
        boxplot[class] = α_vec .+ 1

        # Compute TemperedPareto fit
        agg_df.samples_id .= "samp"
        α, β, ε = fit_samples(agg_df; xmins=xmins)
        # α, β, ε = α[1], β[1], ε[1]
        α, ε = α[1], ε[1]
        fit[class] = MDist.ParetoI(α, ε)
        
        if filter
            mask = counts .> ε
            ranks = ranks[mask]
            counts = counts[mask]
        end

        zipf[class] = (ranks = ranks, counts = counts)
    end

    out = (zipf = zipf, heaps = heaps, boxplot = boxplot, fit = fit)
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
        # β = fit.β
        # push!(β_vec, β)
        ε = fit.ε
        push!(ε_vec, ε)
        α = fit.α
        push!(α_vec, α)
    end
    
    return (α_vec, β_vec, ε_vec)
end

function computevocabsize(df::DataFrame; rng=Random.Xoshiro(42))
    #/ Construct vocabulary and dictionary
    #  note: dictionary is a set for quick comparison
    heapdf = DataFrame(documentsize=Int[], vocabularysize=Int[])
    vocabulary = []
    dictionary = Set()
    #/ In random order, compute the vocabularysize for increasing document sizes
    sample_ids = unique(df[!, :sample_id])
    _order = randperm(rng, length(sample_ids))
    for id in sample_ids[_order]
        idxs = findall(df[!, :sample_id] .== id)
        for word in df[!, :component_id][idxs]
            if !(word in dictionary)
                push!(vocabulary, word)
                push!(dictionary, word)
            end
        end

        _documentsize = isempty(heapdf[!, :documentsize]) ? length(idxs) :
                        last(heapdf[!, :documentsize]) + length(idxs)
        push!(heapdf, [_documentsize, length(vocabulary)])
    end
    return heapdf
end

end # End module