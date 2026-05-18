#= Downsample every standardized dataset and store the results.

By default, each class is downsampled to the smallest sample nreads observed in
that class. Outputs are written to data/downsampled as one JLD2 file per source
dataset. Pass --csv to also write compressed CSV files.
=#

using CSV
using DataFrames
using Dates
using JLD2
using Random

using Meris

const OUTDIR = joinpath(Meris.DATADIR, "downsampled")

const DATASETS = Dict(
    "arxiv" => () -> Meris.arXivLoader.load(),
    "gutenberg" => () -> Meris.GutenbergLoader.load(),
    "rfc" => () -> Meris.RFCLoader.load(),
    "otu" => () -> Meris.OTULoader.load(),
    "finance" => () -> Meris.FinanceLoader.load(),
    "gowalla" => () -> Meris.GowallaLoader.load(),
    "lego" => () -> Meris.LEGOLoader.load(),
    "bci" => () -> Meris.BCITreeLoader.load(),
    "biotime" => () -> Meris.BioTIMELoader.load(),
    "gtex" => () -> Meris.GTExLoader.load(),
)

function parse_args(args)
    selected = sort(collect(keys(DATASETS)))
    write_csv = false
    fail_fast = false
    seed = 123

    for arg in args
        if startswith(arg, "--datasets=")
            selected = split(last(split(arg, "="; limit=2)), ",")
        elseif startswith(arg, "--seed=")
            seed = parse(Int, last(split(arg, "="; limit=2)))
        elseif arg == "--csv"
            write_csv = true
        elseif arg == "--fail-fast"
            fail_fast = true
        elseif arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code julia-code/scripts/cli-scripts/downsample-datasets.jl [options]

            Options:
              --datasets=a,b   Comma-separated subset. Available: $(join(sort(collect(keys(DATASETS))), ", "))
              --seed=N         Random seed for reproducible downsampling. Default: 123
              --csv            Also write .csv.gz files next to the .jld2 files.
              --fail-fast      Stop on the first dataset error.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    unknown = setdiff(selected, keys(DATASETS))
    isempty(unknown) || error("Unknown dataset(s): $(join(unknown, ", "))")

    return (; selected, write_csv, fail_fast, seed)
end

function downsample_summary(df::DataFrame)
    isempty(df.class) && return DataFrame(class=String[], nreads=Int[], nsamples=Int[], ncomponents=Int[], rows=Int[])

    summary = combine(
        groupby(df, :class),
        :nreads => first => :nreads,
        :sample_id => (x -> length(unique(x))) => :nsamples,
        :component_id => (x -> length(unique(x))) => :ncomponents,
        nrow => :rows,
    )
    sort!(summary, :class)

    return summary
end

function save_downsampled(name::AbstractString; write_csv::Bool=false)
    @info "Loading dataset" dataset=name
    df = DATASETS[name]()
    select!(df, :class, :sample_id, :component_id, :counts, :nreads)

    input_summary = combine(
        groupby(df, :class),
        :nreads => minimum => :target_nreads,
        :sample_id => (x -> length(unique(x))) => :nsamples,
        :component_id => (x -> length(unique(x))) => :ncomponents,
        nrow => :rows,
    )
    sort!(input_summary, :class)

    @info "Downsampling dataset" dataset=name classes=nrow(input_summary) rows=nrow(df)
    ds_df = Meris.DataTools.downsample_df(df)
    summary = downsample_summary(ds_df)

    jld2_file = joinpath(OUTDIR, "$(name).jld2")
    dataset = String(name)
    created_at = string(now())
    default_rule = "per-class minimum nreads"
    @save jld2_file dataset ds_df summary input_summary created_at default_rule

    if write_csv
        csv_file = joinpath(OUTDIR, "$(name).csv.gz")
        CSV.write(csv_file, ds_df; compress=true)
    end

    @info "Saved downsampled dataset" dataset=name rows=nrow(ds_df) file=jld2_file

    return (; dataset=name, rows=nrow(ds_df), classes=nrow(summary), file=jld2_file)
end

function main(args=ARGS)
    options = parse_args(args)

    mkpath(OUTDIR)
    Random.seed!(options.seed)

    results = DataFrame(dataset=String[], rows=Int[], classes=Int[], file=String[])
    for name in options.selected
        try
            result = save_downsampled(name; write_csv=options.write_csv)
            push!(results, result)
        catch err
            if options.fail_fast
                rethrow(err)
            end
            @error "Failed to downsample dataset" dataset=name exception=(err, catch_backtrace())
        end
        GC.gc()
    end

    summary_file = joinpath(OUTDIR, "summary.csv")
    CSV.write(summary_file, results)
    @info "Finished downsampling datasets" summary_file=summary_file

    return results
end

main()
