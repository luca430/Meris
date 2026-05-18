#= Downsample standardized datasets grouped as in the plotting scripts.

Groups:
- linguistic: arXiv, Gutenberg, RFC
- microbial: OTU
- social: finance, Gowalla, LEGO
- biology: GTEx, BCI, BioTIME

By default, each class is downsampled to the smallest feasible sample nreads
within that class. Outputs are written to data/downsampled as one JLD2 file per
group. Pass --csv to also write compressed CSV files.
=#

using CSV
using DataFrames
using Dates
using JLD2
using Random

using Meris

const OUTDIR = joinpath(Meris.DATADIR, "downsampled")
const GROUPS = ["linguistic", "microbial", "social", "biology"]

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
    df = Meris.OTULoader.load()
    select!(df, :class, :sample_id, :component_id, :counts, :nreads)
    return df
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

function load_group_df(group::AbstractString)
    group == "linguistic" && return load_linguistic_df()
    group == "microbial" && return load_microbial_df()
    group == "social" && return load_social_df()
    group == "biology" && return load_biology_df()
    error("Unknown group: $group")
end

function parse_args(args)
    selected = copy(GROUPS)
    write_csv = false
    fail_fast = false
    seed = 123

    for arg in args
        if startswith(arg, "--groups=")
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
              julia --project=julia-code julia-code/scripts/cli-scripts/downsample-dataset-groups.jl [options]

            Options:
              --groups=a,b     Comma-separated subset. Available: $(join(GROUPS, ", "))
              --seed=N         Random seed for reproducible downsampling. Default: 123
              --csv            Also write .csv.gz files next to the .jld2 files.
              --fail-fast      Stop on the first group error.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end

    unknown = setdiff(selected, GROUPS)
    isempty(unknown) || error("Unknown group(s): $(join(unknown, ", "))")

    return (; selected, write_csv, fail_fast, seed)
end

function class_summary(df::DataFrame)
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

function input_class_summary(df::DataFrame)
    summary = combine(
        groupby(df, :class),
        :nreads => minimum => :target_nreads,
        :sample_id => (x -> length(unique(x))) => :nsamples,
        :component_id => (x -> length(unique(x))) => :ncomponents,
        nrow => :rows,
    )
    sort!(summary, :class)

    return summary
end

function save_group(group::AbstractString; write_csv::Bool=false)
    @info "Loading group" group
    df = load_group_df(group)
    select!(df, :class, :sample_id, :component_id, :counts, :nreads)

    input_summary = input_class_summary(df)

    @info "Downsampling group" group classes=nrow(input_summary) rows=nrow(df)
    ds_df = Meris.DataTools.downsample_df(df)
    summary = class_summary(ds_df)

    jld2_file = joinpath(OUTDIR, "$(group).jld2")
    created_at = string(now())
    default_rule = "per-class minimum feasible nreads"
    @save jld2_file group ds_df summary input_summary created_at default_rule

    if write_csv
        csv_file = joinpath(OUTDIR, "$(group).csv.gz")
        CSV.write(csv_file, ds_df; compress=true)
    end

    @info "Saved downsampled group" group rows=nrow(ds_df) file=jld2_file

    return (; group, rows=nrow(ds_df), classes=nrow(summary), file=jld2_file)
end

function main(args=ARGS)
    options = parse_args(args)

    mkpath(OUTDIR)
    Random.seed!(options.seed)

    results = DataFrame(group=String[], rows=Int[], classes=Int[], file=String[])
    for group in options.selected
        try
            result = save_group(group; write_csv=options.write_csv)
            push!(results, result)
        catch err
            if options.fail_fast
                rethrow(err)
            end
            @error "Failed to downsample group" group exception=(err, catch_backtrace())
        end
        GC.gc()
    end

    summary_file = joinpath(OUTDIR, "group-summary.csv")
    CSV.write(summary_file, results)
    @info "Finished downsampling groups" summary_file

    return results
end

main()
