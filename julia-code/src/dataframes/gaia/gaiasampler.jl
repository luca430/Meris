#= Module to investigate macro(-ecological) laws in gaia dataset

=#
#/ Start module
module GaiaSampler

#/ Packages
using CSV, Glob, CodecZlib
using DataFrames, DataFramesMeta
using Random, StatsBase

#/ Modules, directories
import Meris.GAIADIR as GAIADIR

#################
### FUNCTIONS ###
"""
Function to import raw data, bin stars and convert to standardized DataFrame
WARNING: this function parses ~40GiB of data so it requires quite some time to conclude
"""
function parse_gaia_raw_data(;
        files = glob("*/sanitized/*.csv.gz", GAIADIR),
        nbins = 200,
        min_b = 0.25,
        max_b = 3.5,
        savefile = false,
        filename = "gaia_data.csv.gz"
    )

    dfs = DataFrame[]
    for f in files
        println(f)
        try
            df = CSV.read(f, DataFrame; select=["BP-RP", "RUWE"])

            filter!(row -> !ismissing(row.RUWE) && row.RUWE isa Real, df)
            filter!(row -> row.RUWE < 1.4, df)
            if nrow(df) < 1
                continue
            end

            df.component_id = star_bp_rp_identifier.(df."BP-RP"; nbins=nbins, min_b=min_b, max_b=max_b)
            df_counts = combine(groupby(df, :component_id), nrow => :counts)
            df_counts.nreads .= sum(df_counts.counts)
            df_counts.sample_id .= extract_sample_id(f)
            filter!(row -> row.component_id > 0, df_counts)
            push!(dfs, df_counts)

        catch e
            @warn "Failed to process file (possibly corrupted): $f" exception=(e, catch_backtrace())
            continue
        end
    end

    final_df = vcat(dfs...)

    if savefile
        path = joinpath(GAIADIR, filename)
        open(GzipCompressorStream, path, "w") do io
            CSV.write(io, final_df)
        end
    end

    return final_df
end

########################
### HELPER FUNCTIONS ###
"Function to extract RA/DEC info from filename"
function extract_sample_id(path)
    m = match(r"ra(\d+-\d+)_dec(-?\d+-?-?\d*)", basename(path))
    return isnothing(m) ? missing : string("ra", m.captures[1], "_dec", m.captures[2])
end

"Function for BP-RP binning"
function star_bp_rp_identifier(bp_rp; nbins=200, min_b=0.25, max_b=3.5)
    Δb = (max_b - min_b) / nbins
    binedges = collect(min_b:Δb:max_b)
    idx = findfirst(binedges .> bp_rp)
    return isnothing(idx) ? length(binedges) + 1 : idx
end

end # module arXivSampler
#/ End module

