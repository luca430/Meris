#= Module to investigate macro(-ecological) laws in genes data from the GTEx project

=#
#/ Start module
module GTExLoader

#/ Packages
using DataFrames, DataFramesMeta
using Random, StatsBase
using CodecZlib, CSV, Glob

#/ Modules, directories
import Meris.GTEXDIR as GTEXDIR

#################
### FUNCTIONS ###
function load(
    ;
    DIR=GTEXDIR * "processed/",
    verbose=false,    
    filterdata    = true,
    minsamples    = 30,
    minreads      = 10^8,
    mincomponents = 100
    )
    files = glob("**/*.long.csv.gz", DIR)
    isempty(files) && error("No *.long.csv.gz files found under $DIR")

    dfs = DataFrame[]

    for f in files
        (verbose) && (println("Loading $f"))
        df = CSV.read(f, DataFrame)
        push!(dfs, df)
    end
    df = vcat(dfs...)
    #~ Rename for consistency
    rename!(df, :tissue => :class, :gene_id => :component_id)
    
    if filterdata
        #~ filter data
        @subset!(df, :nreads .> minreads)
        summarydf = @chain df begin
            @by(
                :class,
                :nsamples = length(:sample_id),
                :ncomponents = length(unique(:component_id))
            )
            @subset(:nsamples .> minsamples, :ncomponents .> mincomponents)
        end
        @subset!(df, :class .∈ Ref(summarydf.class))
    end

    return df
end

### FUNCTIONS TO PROCESS BIG RAW-DATA ###
### ~run `gtex_tree_to_many_csv_gz(GTEXDIR * "processed"; flush_n=200_000)` to process data
tissue_from_path(path::AbstractString) = basename(dirname(path))

function header_cols_gz(path::AbstractString; header_row::Int=3)
    io = GzipDecompressorStream(open(path, "r"))
    for _ in 1:(header_row-1)
        readline(io)
    end
    hdr = chomp(readline(io))
    close(io)

    cols = split(hdr, '\t'; keepempty=true)
    gene_col = String(cols[1])
    sample_cols = String.(cols[2:end])
    return gene_col, sample_cols
end

@inline function split_pad(line::AbstractString, expected_len::Int)
    fields = split(chomp(line), '\t'; keepempty=true)
    if length(fields) < expected_len
        append!(fields, fill("", expected_len - length(fields)))
    end
    return fields
end

function compute_nreads_stream(path::AbstractString; header_row::Int=3)
    _, sample_cols = header_cols_gz(path; header_row=header_row)
    ns = length(sample_cols)
    totals = zeros(Int64, ns)

    io = GzipDecompressorStream(open(path, "r"))
    for _ in 1:header_row
        readline(io)
    end

    expected_len = 1 + ns

    for line in eachline(io)
        fields = split_pad(line, expected_len)
        @inbounds for j in 1:ns
            v = fields[j+1]
            if !isempty(v) && v != "0"
                totals[j] += parse(Int64, v)
            end
        end
    end
    close(io)

    return sample_cols, totals
end

"""
Write long table as gzipped CSV: out_csv_gz should end with .csv.gz
Counts stored as Int64 to avoid overflow.
"""
function write_long_stream_to_csv_gz(path::AbstractString, out_csv_gz::AbstractString;
        header_row::Int=3,
        max_samples::Int=200,
        flush_n::Int=300_000,
        append::Bool=true
    )

    tissue = tissue_from_path(path)

    sample_cols, totals = compute_nreads_stream(path; header_row=header_row)
    ns = minimum([length(sample_cols), max_samples])
    expected_len = 1 + ns

    # open gzip output stream (append not supported for gzip streams in a clean way)
    # so: for each input file, write to its own gz, then you can concatenate later if needed.
    first_write = !append || !isfile(out_csv_gz)
    if append && isfile(out_csv_gz)
        error("Appending to .gz is not recommended. Write one .csv.gz per tissue/file, or overwrite.")
    end

    out_io = GzipCompressorStream(open(out_csv_gz, "w"))
    write(out_io, "gene_id,sample_id,counts,nreads,tissue\n")

    gene_buf   = String[]
    sample_buf = String[]
    counts_buf = Int64[]
    nreads_buf = Int64[]
    tissue_buf = String[]

    function flush!()
        n = length(gene_buf)
        n == 0 && return
        @inbounds for i in 1:n
            write(out_io, gene_buf[i]);     write(out_io, ',')
            write(out_io, sample_buf[i]);   write(out_io, ',')
            write(out_io, string(counts_buf[i])); write(out_io, ',')
            write(out_io, string(nreads_buf[i])); write(out_io, ',')
            write(out_io, tissue_buf[i]);   write(out_io, '\n')
        end
        empty!(gene_buf); empty!(sample_buf); empty!(counts_buf); empty!(nreads_buf); empty!(tissue_buf)
    end

    io = GzipDecompressorStream(open(path, "r"))
    for _ in 1:header_row
        readline(io)
    end

    for line in eachline(io)
        fields = split_pad(line, expected_len)
        gene = replace(fields[1], r"\.\d+$" => "")

        @inbounds for j in 1:ns
            v = fields[j+1]
            if !isempty(v) && v != "0"
                c = parse(Int64, v)
                push!(gene_buf, gene)
                push!(sample_buf, sample_cols[j])
                push!(counts_buf, c)
                push!(nreads_buf, totals[j])
                push!(tissue_buf, tissue)
            end
        end

        if length(gene_buf) >= flush_n
            flush!()
        end
    end

    close(io)
    flush!()
    close(out_io)

    return out_csv_gz
end

function gtex_tree_to_many_csv_gz(
    root::AbstractString;
    header_row::Int=3,
    max_samples::Int=200,
    flush_n::Int=200_000,
    verbose=false                             
    )

    files = glob("**/*.nonzero_cols.gz", root)
    isempty(files) && error("No *.nonzero_cols.gz under $root")

    outs = String[]
    for f in files
        (verbose) && (println("[file] ", f))
        #~ Build output path in SAME directory
        out = replace(f, ".nonzero_cols.gz" => ".long.csv.gz")
        #~ Write to DataFrame
        push!(outs,
            write_long_stream_to_csv_gz(
                f,
                out;
                header_row=header_row,
                max_samples=max_samples,
                flush_n=flush_n,
                append=false
            )
        )
        (verbose) && (println("  -> ", out))
    end
    return outs
end


end # module GTExLoader
#/ End module

