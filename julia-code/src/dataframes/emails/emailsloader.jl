#= Module to load email dataset

=#
#/ Start module
module EmailsLoader

#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase

#/ Modules, directories
import Meris.EMAILDIR as EMAILDIR

#################
### FUNCTIONS ###
"""
Build counts DataFrame: sample_id (Date line raw string), component_id, counts
"""
function load(rootdir::AbstractString)
    sample_ids = String[]
    component_ids = String[]

    for (dir, _, files) in walkdir(rootdir)
        # only directories that are inside "_sent_mail"
        if occursin(string("/", "_sent_mail"), dir)
            for f in files
                path = joinpath(dir, f)
                comp = splitpath(path)[end-2]  # expects rootdir ends with "maildir"
                date = extract_date_line(path)[end-7:end]

                if comp !== missing && date !== missing
                    push!(component_ids, String(comp))
                    push!(sample_ids, String(date))
                end
            end
        end
    end

    df = DataFrame(sample_id = sample_ids, component_id = component_ids)
    gdf = groupby(df, [:sample_id, :component_id])
    out = combine(gdf, nrow => :counts)
    out.class .= "ENRON"
    return @chain out begin
        @groupby(:sample_id)
        @combine(:class, :sample_id, :component_id, :counts, :nreads=sum(:counts))
    end
end

### HELPER ###
"""
Extract "day month year" from Date header using only string operations.
Example output: "14 May 2001"
"""
function extract_date_line(path::AbstractString)
    open(path, "r") do io
        for line in eachline(io)
            if startswith(line, "Date:")
                raw = strip(line[length("Date:")+1:end])

                # remove weekday if present ("Mon, ")
                raw = replace(raw, r"^[A-Za-z]+,\s*" => "")

                parts = split(raw)

                if length(parts) >= 3
                    return string(parts[1], " ", parts[2], " ", parts[3])
                else
                    return missing
                end
            end
        end
    end
end

end # module GowallaLoader
#/ End module

