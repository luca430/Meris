#= Module to sample from the OTU dataset
   OTU data can be obtained from https://github.com/jacopogrilli/lawsdiv,
   which has an `.RData` file that contains all data that is needed for analysis.
=#
#/ Start module
module OTULoader

#/ Packages
using CSV
using DataFrames, DataFramesMeta
using Random
using RData
using StatsBase

using Meris

#/ Modules, directories
import Meris.OTUDIR as OTUDIR

const ENV_MAP = Dict(
    " seawater.MGYS00002437" => "SEA",
    "ORAL.SRP056641" => "ORAL1",
    "GUT.SRP056641" => "GUT1",
    "VAGINAL.SRP056641" => "VAGINAL",
    "Environmental Terrestrial Soil.SRP052295" => "SOIL",
    "oralcavity.ERP021896" => "ORAL2",
    "feces.ERP021896" => "FECES",
    "skin.ERP021896" => "SKIN",
    "Environmental Aquatic Marine Hydrothermal vents.ERP017354" => "AQUA1",
    "GUT.ERP015450" => "GUT2",
    "GUT.ERP013827" => "GUT3",
    "River.ERP012927" => "RIVER",
    "Lake.ERP012927" => "LAKE",
    "Environmental Aquatic Marine.ERP009703" => "AQUA2",
    "activatedsludge.ERP009143" => "SLUDGE",
    "Glacier.ERP017997" => "GLACIER"
)

#################
### FUNCTIONS ###
"""
    load_rdata

Load RData file into a DataFrame
"""
function load_rdata(; rdatafilename = ROTUDIR*"crosssecdata.RData")
    df = RData.load(rdatafilename)["datatax"]
    df = @transform(df, :classification = String.(:classification))
    return df
end

"Convert raw DataFrame to the standardized version"
function standardized(df)
    df = @chain df begin
        @transform(:class = :classification .* "." .* :project_id)
        @rename(:counts = :count, :component_id = :otu_id)
        @select(:class, :component_id, :counts, :nreads, :run_id)
        @rename(:sample_id = :run_id)
    end
    df.class = get.(Ref(ENV_MAP), df.class, df.class)
    return df
end

"Load all EBI Metagenomics OTU data into a single DataFrame"
function load(
    ;
    datafilename = OTUDIR * "/raw-data/crosssecdata.RData",    
    filterdata    = true,
    minreads::Int=10_000,
    mincomponents::Int=500,
    minsamplecomponents::Int=200,
    minsamples::Int=30,   
    )
    df = load_rdata(; rdatafilename = datafilename)
    df = standardized(df)
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

end # module OTUSampler
#/ End module
