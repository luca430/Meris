module SAD

using Meris
using DataFrames, DataFramesMeta, StatsBase
using JLD2, Random
import Meris.MDistributions as MDist

function compute(
        df;
        pareto=:ParetoI,
        nbins=30,
        save=true,
        filename="SAD.jld2"
    )

    Random.seed!(1234)
    classes = unique(df.class)

    heaps_d = Dict()
    pl_d = Dict()
    ple_d = Dict()
    cad_d = Dict()
    for class in classes
        println(class)
        sdf = df[df.class .== class, :]
    
        # Compute Heaps' scaling
        hdf = heaps(sdf)
        heaps_d[class] = (N = hdf.documentsize, V = hdf.vocabularysize)
    
        # Compute power law exponents
        α_vec, ε_vec = get_params(sdf; pareto=pareto)
        pl_d[class] = (α = α_vec, ε = ε_vec)
    
        # Compute CAD
        agg_df = aggregate_samples(sdf)
        agg_df.counts .= agg_df.counts ./ agg_df.nreads
        x = log10.(agg_df.counts[agg_df.counts .> mean(ε_vec)])
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
        @combine(:sample_id, counts = sum(:counts), nreads = sum(:nreads))
    end
end

function get_params(df; pareto=:ParetoI)
    a_idx, e_idx = pareto_idx(pareto)
    a_vec, e_vec = Float64[], Float64[]
    for row in eachrow(df)
        params = row[pareto]
        a = params[a_idx]
        e = params[e_idx]
        if pareto == :GeneralizedPareto
            a = 1/a
        end
        push!(a_vec, a)
        push!(e_vec, e)
    end

    return a_vec, e_vec
end

function pareto_idx(pareto)
    if pareto == :ParetoI
        return 1, 2
    elseif pareto == :TemperedPareto
        return 1, 3
    elseif pareto == :GeneralizedPareto
        return 3, 1
    elseif pareto == :ParetoIV
        return 1, 4
    end
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