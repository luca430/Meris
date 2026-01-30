#= Module to investigate macro(-ecological) laws in gaia dataset

=#
#/ Start module
module GaiaLoader

#/ Packages
using CSV, Glob, CodecZlib, CategoricalArrays
using DataFrames, DataFramesMeta, StatsBase

#/ Modules, directories
import Meris.GAIADIR as GAIADIR

const COLOR_BINS = 0.5:0.05:4.5
const DIST_BINS  = 0:100:20000

#################
### FUNCTIONS ###
function load(; file = GAIADIR * "processed/gaia_binned_counts.csv.gz")
    df = CSV.read(file, DataFrame)
    df.class .= "GAIA"
    df.component_id .= df.color_bin .* df.dist_bin

    sdf = @chain df begin
        @groupby(:sample)
        @combine(
            :class,
            :component_id,
            :count,
            :nreads = sum(:count)
        )
    end

    rename!(sdf, :count => :counts, :sample => :sample_id)
    return sdf
end

### DATA CONVERSION ###
"""
Functions to import raw data, bin stars and convert to standardized DataFrame
WARNING: this function parses ~40GiB of data so it requires quite some time to conclude
"""
function sample_key(path::AbstractString)
    m = match(r"_ra\d+(?:\.\d+)?-\d+(?:\.\d+)?_dec-?\d+(?:\.\d+)?-?-?\d+(?:\.\d+)?_", path)
    m === nothing && error("Cannot parse sample from filename: $path")
    strip(m.match, '_')
end

function safe_read_bprp_rgeo(path::AbstractString)
    empty = DataFrame("BP-RP" => Float64[], "rgeo" => Float64[])

    df = try
        CSV.read(path, DataFrame;
            select = ["BP-RP", "rgeo"],
            types  = Dict("BP-RP" => String, "rgeo" => String),
            strict = false,
            silencewarnings = true,
            validate = false
        )
    catch err
        @warn "CSV.read failed, skipping file" path exception=(err, catch_backtrace())
        return empty
    end

    # if CSV saw no header / no cols, skip
    if !("BP-RP" in names(df)) || !("rgeo" in names(df))
        @warn "Missing required columns, skipping file" path
        return empty
    end

    df."BP-RP" = [x === missing ? missing : tryparse(Float64, String(x)) for x in df."BP-RP"]
    df.rgeo    = [x === missing ? missing : tryparse(Float64, String(x)) for x in df.rgeo]
    
    dropmissing!(df)
    return df
end

function counts_table_for_file(path::AbstractString;
        color_bins=COLOR_BINS, dist_bins=DIST_BINS)

    df = safe_read_bprp_rgeo(path)
    select!(df, "BP-RP", "rgeo")
    dropmissing!(df)
    df = df[(df."BP-RP" .> minimum(COLOR_BINS)) .&& (df."BP-RP" .< maximum(COLOR_BINS)),:]
    df = df[(df.rgeo .> minimum(DIST_BINS)) .&& (df.rgeo .< maximum(DIST_BINS)),:]

    df.color_bin = cut(df."BP-RP", color_bins)
    df.dist_bin  = cut(df.rgeo, dist_bins)
    filter!(r -> !ismissing(r.color_bin) && !ismissing(r.dist_bin), df)

    ct = combine(groupby(df, [:color_bin, :dist_bin]), nrow => :count)
    ct.sample .= sample_key(path)

    @info "Parsed" file=basename(path) n=nrow(df)
    return ct
end

function counts_all_samples(dir::AbstractString; pattern="gaia_stars_*_cleaned.csv.gz")
    paths = sort(glob(pattern, dir))
    isempty(paths) && error("No files matched in $dir with pattern $pattern")
    df = vcat((counts_table_for_file(p) for p in paths)...)

    return filter!(r -> !ismissing(r.color_bin) && !ismissing(r.dist_bin), df)
end

end # module arXivSampler
#/ End module

