#= Module to load gutenberg dataset

=#
#/ Start module
module GutenbergLoader

#/ Packages
using Glob
using CSV, DataFrames, DataFramesMeta
using Random, StatsBase

#/ Modules, directories
import ..DataTools: filterdata
import Meris.GUTENBERGDIR as GUTENBERGDIR

#################
### FUNCTIONS ###
function load(
    ;
    root=GUTENBERGDIR * "raw-data",
    marker=r"\*\*\*.*\*\*\*",
    minsamples    = 30,
    minreads      = 100_000,
    mincomponents = 100,
    applyfilter   = true,
    reorder       = true,
    top           = nothing,
    )
    df = DataFrame(
        class=String[],
        component_id=String[],
        sample_id=String[],
        counts=Int[],
        nreads=Int[]
    )
    #~ Load data
    for (dirpath, _dirs, files) in walkdir(root)
        for f in files
            endswith(lowercase(f), ".txt") || continue
            path = joinpath(dirpath, f)

            cls = class_from_path(path, root)
            sid = sample_id_from_path(path)

            text = load_text(path)
            text === nothing && continue   # skip this file
            
            tokens = tokenize(text)
            cnt_map = countmap(tokens)
            nreads = length(tokens)
            
            for (tok, c) in cnt_map
                push!(df, (cls, tok, sid, c, nreads))
            end
        end
    end

    if applyfilter
        #~ filter data
        df = filterdata(
            df; minsamples=minsamples, minreads=minreads, mincomponents=mincomponents,
            reorder=reorder, top=top
        )
    end

    return df
end

##############
### HELPER ###
function load_text(path::AbstractString; marker=r"\*\*\*.*\*\*\*")
    lines = readlines(path)
    idx = findfirst(l -> occursin(marker, l), lines)

    idx === nothing && return nothing   # skip file

    join(lines[(idx+1):end], "\n")
end

function tokenize(text::String)
    text = lowercase(text)
    text = replace(text, r"[^\p{L}\s]+" => " ")
    split(text)
end

# --- helpers ---
"Class is the first folder under root, e.g. root/en/..../file.txt -> 'en'"
function class_from_path(path::AbstractString, root::AbstractString)
    rel = relpath(path, root)
    parts = splitpath(rel)
    isempty(parts) && error("Could not get class from path: $path")
    return parts[1]
end

"Sample id from filename: '077056_The_dream_detective.txt' -> '077056' (fallback to stem)"
function sample_id_from_path(path::AbstractString)
    base = basename(path)                    # e.g. 077056_The_dream_detective.txt
    stem = splitext(base)[1]                 # e.g. 077056_The_dream_detective
    m = match(r"^\d+", stem)                 # take leading digits if present
    return m === nothing ? stem : m.match
end

end # module GutenbergSampler
#/ End module

