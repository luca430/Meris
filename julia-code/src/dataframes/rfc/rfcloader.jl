#= Module to investigate macro(-ecological) laws in RFC document =#
#/ Start module
module RFCLoader

#/ Packages
using CSV, DataFrames, DataFramesMeta
using Random, StatsBase

using Meris

#/ Modules, directories
import Meris.RFCDIR as RFCDIR

#################
### FUNCTIONS ###
"Collect all tokens of all RFCs in a single DataFrame for post-processing"
function load(
    ;
    DIR=RFCDIR * "raw-data/",
    maxfiles      = 5000,      #~ Max. no of parsed files
    maxrows       = 2_000_000, #~ Max. no of rows allowed in DataFrame
    filterdata    = true,
    minreads::Int=10_000,
    mincomponents::Int=1_000,
    minsamplecomponents::Int=500,
    minsamples::Int=30,           
)
    #~ Allocate DataFrame
    df = DataFrame(
        sample_id=String[], component_id=String[], counts=Int[], nreads=Int[], class=String[]
    )
    nfiles = 0
    #~ Loop through all files in DIR, and extract tokens
    for FILE in readdir(DIR)
        #~ Check if filename is `rfc[0-9].txt`
        if occursin(r"rfc[0-9]+\.txt$", basename(FILE))
            tokens = load_rfc(; FILENAME=FILE, DIR=DIR)
            if length(tokens) > minsamplecomponents
                nfiles += 1
                #/ Perform a simple countmap
                cm = countmap(tokens)
                nreads = length(tokens)
                #/ Put tokens in DataFrame
                for (token, count) in pairs(cm)
                    push!(df, [splitext(FILE)[begin], token, count, nreads, "rfc"])
                end
            end
        end
        #~ Stop if max. files have been tokenized or if maximum no. of rows is reached
        (has_reached(nfiles, maxfiles) || has_reached(nrow(df), maxrows)) && (break)
    end

    if filterdata
        #~ filter data
        df = Meris.DataTools.df_filter(
            df;
            minreads=minreads,
            mincomponents=mincomponents,
            minsamplecomponents=minsamplecomponents,
            minsamples=minsamples
        )
    end
    
    return df
end

########################
### HELPER FUNCTIONS ###
"Load RFC document"
function load_rfc(;
    FILENAME="rfc1.txt",
    DIR=RFCDIR * "raw-data/"
)
    txt = read(DIR * FILENAME, String)
    tokens = collect(eachmatch(r"[A-Za-z]+", txt)) .|> x -> lowercase(x.match)
    return tokens
end

function has_reached(x, y)
    return !isnothing(y) && (x > y)
end

end # module RFCLoader
#/ End module
