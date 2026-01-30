#= Module to investigate macro(-ecological) laws in gaia dataset

=#
#/ Start module
module GaiaLoader

#/ Packages
using CSV, Glob, CodecZlib, CategoricalArrays
using DataFrames, DataFramesMeta, StatsBase

#/ Modules, directories
import Meris.GAIADIR as GAIADIR

const COLOR_BINS = 0.5:0.05:4.0
const DIST_BINS  = 0:200:30000

#################
### FUNCTIONS ###
function load(; file = GAIADIR * "processed/clean_F-G-K-M_nbin_500_bmin_0.25_bmax_3.5.csv.gz")
    df = CSV.read(file, DataFrame)
    df.class .= "gaia"
    return df
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
    df = CSV.read(path, DataFrame;
        select = ["BP-RP", "rgeo"],          # only read what you need
        types  = Dict("BP-RP" => String, "rgeo" => String),  # read as strings first
        strict = false,
        silencewarnings = true,
    )

    # robust conversion: non-parsable -> missing
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

    df.color_bin = cut(df."BP-RP", color_bins)
    df.dist_bin  = cut(df.rgeo, dist_bins)
    filter!(r -> !ismissing(r.color_bin) && !ismissing(r.dist_bin), df)

    ct = combine(groupby(df, [:color_bin, :dist_bin]), nrow => :count)
    ct.sample .= sample_key(path)
    return ct
end

function counts_all_samples(dir::AbstractString; pattern="gaia_stars_*_cleaned.csv.gz")
    paths = sort(glob(pattern, dir))
    isempty(paths) && error("No files matched in $dir with pattern $pattern")

    return vcat((counts_table_for_file(p) for p in paths)...)
end

end # module arXivSampler
#/ End module

