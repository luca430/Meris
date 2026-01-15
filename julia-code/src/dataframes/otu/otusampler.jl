#= Module to sample from the OTU dataset
   OTU data can be obtained from https://github.com/jacopogrilli/lawsdiv, which has an `.RData`
   file that contains the information that is needed.
=#
#/ Start module
module OTUSampler

#/ Packages
using CSV
using DataFrames, DataFramesMeta
using Random
using RData
using StatsBase

#/ Modules, directories
import Meris.OTUDIR as OTUDIR
const ROTUDIR = OTUDIR * "RData/"
const CSVOTUDIR = OTUDIR * "csv/"

#################
### FUNCTIONS ###
function samplevocabsize(
    edf::DataFrame,  #~ DataFrame containing the OTUs of a specific environment
    Nv::Vector{Int};
    n::Int = 32,
    rng = Random.Xoshiro(42*n),
    filterdf = true
)
    (filterdf) && (edf = filter_data(edf))
	  #~ Create "bag of OTUs"
    bagofotus = @select(edf, :otu_id, :counts)
    w = Weights(bagofotus.counts)
    #~ For each N in Nv, sample N OTUs n times and compute the vocabulary size V(N)
    V = zeros(Int, length(Nv), n)
    for i in eachindex(Nv), k in 1:n
        s = StatsBase.sample(rng, bagofotus.otu_id, w, Nv[i], replace=false)
        V[i,k] = length(unique(s))
    end
    return V
end


########################
### HELPER FUNCTIONS ###
"""
    load_rdata

Load RData file into a DataFrame
"""
function load_rdata(; rdatafilename = ROTUDIR*"crosssecdata.RData")
    @info "Loading raw RData..."
    df = RData.load(rdatafilename)["datatax"]
    df = @transform(df, :classification = String.(:classification))
    return df
end

"""
    load_data

Load DataFrame from specific environment
"""
function load_data(env; filename = CSVOTUDIR * "rawotudata_$(env).csv")
	  return CSV.read(filename, DataFrame, delim=", ")
end

"""
    split_data

Split data into data for distinct environments, such that for each environment we
have a separate database that can be analyzed
"""
function split_data(df::DataFrame; dry=false)
    @info "Splitting raw data..."
    #/ 1. Create unique ID for project and classification
    df = @transform(df, :projectclassification = :project_id .* :classification)
    #~ Load short-hand names for the environments
    #  these are defined in CSVDATAPATH/environmentnames.csv
    environmentnamesdf = CSV.read(
        CSVOTUDIR*"environmentnames.csv", DataFrame, delim=", ",
        types=Dict(:projectclassification => String, :environmentname => String)
    )
    #~ Replace the :project_id and :classification columns by a single column
    #  this makes it easier down the line, and also provides the :environmentname column
    df = @chain begin
        innerjoin(df, environmentnamesdf, on=:projectclassification)
        select(Not([:project_id, :classification, :projectclassification]))
    end

    #/ 2. For each environmentname, select the subset of data from that environment and
    #     store it seperately for further analysis
    for environment in environmentnamesdf.environmentname
        edf = @subset(df, :environmentname .== environment)
        #~ Convert some default types to standard or more appropriate ones
        transform!(edf, :count => ByRow(Int) => :counts)
        transform!(edf, :nreads => ByRow(Int) => :nreads)
        edf = @select(edf, Not([:count]))
        #~ Save dataframe
        if !dry 
            filename = CSVOTUDIR * "rawotudata_$(environment).csv"
            CSV.write(filename, edf, delim=", ")
        end        
    end
    return environmentnamesdf
end

"""
    filter_data

Filter DataFrame to include only
- environments with sufficient samples
- experiments with sufficient reads           [total no. of hits for specific species]
- experiments with sufficient vocabulary size [no. of distinct species]
- species with sufficient counts              [hits of a specific species]

Species are determined on the OTU level
"""
function filter_data(df::DataFrame;
    minsamples = 30,
    minreads = 10_000,
    mincounts = 30,
    minspecies = 50,
    remove_runs = ["ERR1104477", "ERR1101508", "SRR2240575"] # bad runs filtered by Grilli
)
    #/ Check the total number of samples
    #~ if not enough samples, do nothing (i.e., skip it)
    nsamples = length(unique(df[!,:sample_id]))
    if nsamples < minsamples
        @info "not enough samples, skipping [env: $(first(df[!,:environmentname]))]"
        return nothing
    end

    #/ Filter entries and create filtered dataframe
    fdf = @chain df begin
        @rsubset(!in(:run_id, remove_runs))
        @subset(:nreads .> minreads)
        @subset(:counts .> mincounts)
    end
    sdf = @by(fdf, :sample_id, :nspecies = length(unique(:otu_id)))
    fdf = @chain begin
        innerjoin(fdf, sdf, on=:sample_id)
        @subset(:nspecies .> minspecies)
    end
    #/ Return filtered dataframe only when non-empty
    (nrow(fdf) > 0) && (return fdf)
    #~ Otherwise, return nothing
    @info "not enough reads or counts [env: $(first(df[!,:environmentname]))]"
    return nothing
end

end # module OTUSampler
#/ End module
