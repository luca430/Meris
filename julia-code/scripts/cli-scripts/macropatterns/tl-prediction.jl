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
DOWNSAMPLED_DIR = Meris.DATADIR * "downsampled/"

#~ Specify variables
save = true
categories = ["linguistic", "microbial", "social", "biology"]
occ = 0.05
frac = 0.5
nreps = 100
seed = 123
min_reps = 10
nbins = 60
min_bin_components = 20

function parse_args!(args)
    global categories, seed, occ, frac, nbins, min_bin_components, save, DOWNSAMPLED_DIR, RESULTDIR

    for arg in args
        if startswith(arg, "--categories=")
            categories = split(last(split(arg, "="; limit=2)), ",")
        elseif startswith(arg, "--downsampled-dir=")
            DOWNSAMPLED_DIR = last(split(arg, "="; limit=2))
        elseif startswith(arg, "--result-dir=")
            RESULTDIR = last(split(arg, "="; limit=2))
        elseif startswith(arg, "--seed=")
            seed = parse(Int, last(split(arg, "="; limit=2)))
        elseif startswith(arg, "--occ=")
            occ = parse(Float64, last(split(arg, "="; limit=2)))
        elseif startswith(arg, "--frac=")
            frac = parse(Float64, last(split(arg, "="; limit=2)))
        elseif startswith(arg, "--nbins=")
            nbins = parse(Int, last(split(arg, "="; limit=2)))
        elseif startswith(arg, "--min-bin-components=")
            min_bin_components = parse(Int, last(split(arg, "="; limit=2)))
        elseif arg == "--no-save"
            save = false
        elseif arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code julia-code/scripts/cli-scripts/macropatterns/tl-prediction.jl [options]

            Inputs:
              Reads grouped downsampled datasets from $(DOWNSAMPLED_DIR)<category>.jld2.
              Create them with scripts/cli-scripts/downsample-dataset-groups.jl.

            Options:
              --categories=a,b           Comma-separated subset. Default: linguistic,microbial,social,biology
              --downsampled-dir=DIR      Directory with grouped downsampled .jld2 files.
                                         Default: $(DOWNSAMPLED_DIR)
              --result-dir=DIR           Directory for TL-prediction outputs.
                                         Default: $(RESULTDIR)
              --seed=N                   Random seed. Default: 123
              --occ=X                    Occupancy filter. Default: 0.05
              --frac=X                   Train fraction. Default: 0.5
              --nbins=N                  Number of coefficient bins. Default: 60
              --min-bin-components=N     Minimum components for selected bins. Default: 20
              --no-save                  Run without writing result files.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    DOWNSAMPLED_DIR = abspath(DOWNSAMPLED_DIR)
    RESULTDIR = abspath(RESULTDIR)
    mkpath(RESULTDIR)
end

function fit_component_bins(test_dfs::Dict, train_summaries::Dict, fit_summaries::Dict)
    dfs = DataFrame[]

    for class in sort(collect(keys(test_dfs)); by=string)
        test_summary = select(
            TLPrediction.summary_df(test_dfs[class]),
            :component_id,
            :mean => :mean_B,
            :var => :var_B
        )
        train_summary = select(
            train_summaries[class],
            :component_id,
            :mean => :mean_A,
            :var => :var_A,
            :coeff,
            :coeff_bin,
            :bin_left_log,
            :bin_right_log,
            :bin_center_log,
            :bin_left,
            :bin_right,
            :bin_center,
            :ncomponents
        )

        component_bins = innerjoin(test_summary, train_summary, on=:component_id)
        component_bins = innerjoin(component_bins, fit_summaries[class], on=:coeff_bin)
        component_bins.class = fill(string(class), nrow(component_bins))
        component_bins.coeff_order = floor.(Int, component_bins.bin_center_log)
        push!(dfs, component_bins)
    end

    return vcat(dfs...)
end

function selected_order_bins(component_bins::DataFrame; min_components::Int=10)
    bin_counts = combine(
        groupby(component_bins, [:class, :coeff_order, :coeff_bin, :bin_center_log, :bin_center, :C_fit_A, :C_fit_B]),
        nrow => :ncomponents
    )
    filter!(:ncomponents => >=(min_components), bin_counts)
    filter!(:C_fit_A => c -> isfinite(c) && c > 0, bin_counts)
    filter!(:C_fit_B => c -> isfinite(c) && c > 0, bin_counts)

    isempty(bin_counts.coeff_order) && return bin_counts

    selected = combine(groupby(bin_counts, [:class, :coeff_order])) do df
        df[argmax(df.ncomponents), :]
    end
    sort!(selected, [:class, :coeff_order])

    return selected
end

function load_category_df(category::AbstractString)
    filename = joinpath(DOWNSAMPLED_DIR, "$(category).jld2")
    isfile(filename) || error(
        "Missing downsampled dataset: $filename. " *
        "Run julia --project=julia-code julia-code/scripts/cli-scripts/downsample-dataset-groups.jl first."
    )

    jldb = JLD2.load(filename)
    haskey(jldb, "ds_df") || error("Downsampled file does not contain ds_df: $filename")

    df = jldb["ds_df"]
    select!(df, :class, :sample_id, :component_id, :counts, :nreads)

    return df
end

function prepare_category(category::AbstractString)
    @info "Loading downsampled $category data..."
    df = load_category_df(category)

    @info "Preparing $category TL prediction data..."
    class_dfs = TLPrediction.divide_in_class(df; downsample=false, occ=occ)
    train_dfs, test_dfs = TLPrediction.partition_per_class(
        class_dfs;
        rng=MersenneTwister(seed),
        frac=frac
    )

    train_summaries = TLPrediction.train_summary_direct(
        train_dfs;
        nbins=nbins
    )
    train_fit_summaries = TLPrediction.train_fit_summary(train_summaries)
    test_fit_summaries = TLPrediction.test_fit_summary(test_dfs, train_summaries)
    fit_summaries = TLPrediction.fit_comparison_summary(train_fit_summaries, test_fit_summaries)
    plot_summaries = fit_summaries

    component_bins = fit_component_bins(test_dfs, train_summaries, fit_summaries)
    selected_bins = selected_order_bins(component_bins; min_components=min_bin_components)

    if save
        filename = joinpath(RESULTDIR, "$(category).jld2")
        downsampled_file = joinpath(DOWNSAMPLED_DIR, "$(category).jld2")
        @save filename train_summaries train_fit_summaries test_fit_summaries fit_summaries plot_summaries component_bins selected_bins min_bin_components occ frac nreps seed min_reps nbins downsampled_file
    end

    df = class_dfs = train_dfs = test_dfs = nothing
    GC.gc()

    return nothing
end

parse_args!(ARGS)

for category in categories
    prepare_category(category)
end
