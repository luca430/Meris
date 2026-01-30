#= Module to sample from the OTU dataset
   OTU data can be obtained from https://github.com/jacopogrilli/lawsdiv, which has an `.RData`
   file that contains the information that is needed.
=#
#/ Start module
module OTULoader

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
    @info "Loading raw RData..."
    df = RData.load(rdatafilename)["datatax"]
    df = @transform(df, :classification = String.(:classification))
    return df
end

"Convert raw DataFrame to the standardized version"
function standardized(df)
    df = @chain df begin
        @transform(:environment = :classification .* "." .* :project_id)
        @rename(:counts = :count, :component_id = :otu_id)
        @select(:environment, :component_id, :counts, :nreads, :run_id)
        @rename(:sample_id = :run_id)
    end
    df.environment = get.(Ref(ENV_MAP), df.environment, df.environment)
    return df
end

function load(; datafilename = ROTUDIR*"crosssecdata.RData")
    df = load_rdata(; rdatafilename = datafilename)
    return standardized(df)
end

end # module OTUSampler
#/ End module
