#= Prepare Taylor's-law prediction summaries for the Figure 2 dataset groups. =#

#/ Packages
using DataFrames
using JLD2
using Random

#/ Modules
using Meris
include("../../module-scripts/macropatterns/TL_prediction.jl")
using .TLPrediction

RESULTDIR = Meris.DATADIR * "macro/tl-prediction/"
mkpath(RESULTDIR)

#~ Specify variables
save = true
# categories = ["linguistic", "microbial", "social", "biology"]
categories = ["social", "biology"]
occ = 0.05
frac = 0.5
nreps = 100
seed = 123
min_reps = 10
nbins = 100
min_bin_components = 10

function test_component_bins(test_dfs::Dict, train_summaries::Dict, test_summaries::Dict)
    dfs = DataFrame[]

    for class in sort(collect(keys(test_dfs)); by=string)
        test_summary = select(
            TLPrediction.summary_df(test_dfs[class]),
            [:component_id, :mean, :var]
        )
        component_bins = innerjoin(
            test_summary,
            train_summaries[class],
            on=:component_id
        )
        component_bins = innerjoin(
            component_bins,
            select(test_summaries[class], :coeff_bin, :C_est, :C_est_err, :C_fit, :C_fit_err),
            on=:coeff_bin
        )
        component_bins.class = fill(string(class), nrow(component_bins))
        component_bins.coeff_order = floor.(Int, component_bins.bin_center_log)
        push!(dfs, component_bins)
    end

    return vcat(dfs...)
end

function selected_order_bins(component_bins::DataFrame; min_components::Int=10)
    bin_counts = combine(
        groupby(component_bins, [:class, :coeff_order, :coeff_bin, :bin_center_log, :bin_center, :C_est, :C_fit]),
        nrow => :ncomponents
    )
    filter!(:ncomponents => >=(min_components), bin_counts)
    filter!(:C_fit => c -> isfinite(c) && c > 0, bin_counts)
    filter!(:C_est => c -> isfinite(c) && c > 0, bin_counts)

    isempty(bin_counts.coeff_order) && return bin_counts

    selected = combine(groupby(bin_counts, [:class, :coeff_order])) do df
        df[argmax(df.ncomponents), :]
    end
    sort!(selected, [:class, :coeff_order])

    return selected
end

function load_linguistic_df()
    df_arxiv = Meris.arXivLoader.load()
    df_arxiv.class .= "arx-" .* uppercase.(df_arxiv.class)
    select!(df_arxiv, :class, :sample_id, :component_id, :counts, :nreads)

    df_gut = Meris.GutenbergLoader.load()
    df_gut.class = "guten-" .* uppercase.(df_gut.class)
    select!(df_gut, :class, :sample_id, :component_id, :counts, :nreads)

    df_rfc = Meris.RFCLoader.load()
    df_rfc.class .= uppercase.(df_rfc.class)
    select!(df_rfc, :class, :sample_id, :component_id, :counts, :nreads)

    return vcat(df_arxiv, df_gut, df_rfc)
end

function load_microbial_df()
    df_otu = Meris.OTULoader.load()
    select!(df_otu, :class, :sample_id, :component_id, :counts, :nreads)
    return df_otu
end

function load_social_df()
    df_fin = Meris.FinanceLoader.load()
    df_fin = df_fin[endswith.(df_fin.class, "-daily"), :]
    df_fin.class = replace.(df_fin.class, "-daily" => "")
    df_fin.class .= "stock-" .* uppercase.(df_fin.class)
    select!(df_fin, :class, :sample_id, :component_id, :counts, :nreads)

    df_gow = Meris.GowallaLoader.load()
    df_gow.class .= "CHECK-IN"
    select!(df_gow, :class, :sample_id, :component_id, :counts, :nreads)

    df_lego = Meris.LEGOLoader.load()
    df_lego.class .= "LEGO"
    select!(df_lego, :class, :sample_id, :component_id, :counts, :nreads)

    return vcat(df_fin, df_gow, df_lego)
end

function load_biology_df()
    df_bci = Meris.BCITreeLoader.load()
    df_bci.class .= "eco-BCI"
    select!(df_bci, :class, :sample_id, :component_id, :counts, :nreads)

    df_bio = Meris.BioTIMELoader.load()
    df_bio.class .= "eco-BT" .* string.(df_bio.class)
    select!(df_bio, :class, :sample_id, :component_id, :counts, :nreads)

    df_gtex = Meris.GTExLoader.load()
    df_gtex.class .= "gen-" .* string.(df_gtex.class)
    select!(df_gtex, :class, :sample_id, :component_id, :counts, :nreads)

    return vcat(df_gtex, df_bci, df_bio)
end

function load_category_df(category::AbstractString)
    category == "linguistic" && return load_linguistic_df()
    category == "microbial" && return load_microbial_df()
    category == "social" && return load_social_df()
    category == "biology" && return load_biology_df()
    error("Unknown category: $category")
end

function prepare_category(category::AbstractString)
    @info "Loading $category data..."
    df = load_category_df(category)

    @info "Preparing $category TL prediction data..."
    class_dfs = TLPrediction.divide_in_class(df; downsample=true, occ=occ)
    train_dfs, test_dfs = TLPrediction.partition_per_class(
        class_dfs;
        rng=MersenneTwister(seed),
        frac=frac
    )

    train_summaries = TLPrediction.train_summary(
        train_dfs;
        nreps=nreps,
        seed=seed,
        frac=frac,
        min_reps=min_reps,
        nbins=nbins
    )
    test_summaries = TLPrediction.test_summary(test_dfs, train_summaries)

    component_bins = test_component_bins(test_dfs, train_summaries, test_summaries)
    selected_bins = selected_order_bins(component_bins; min_components=min_bin_components)

    if save
        filename = RESULTDIR * "$(category).jld2"
        @save filename train_summaries test_summaries component_bins selected_bins min_bin_components occ frac nreps seed min_reps nbins
    end

    df = class_dfs = train_dfs = test_dfs = nothing
    GC.gc()

    return nothing
end

for category in categories
    prepare_category(category)
end
