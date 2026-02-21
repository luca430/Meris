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
    @info "Loading raw RData..."
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
    minsamples    = 30,
    minreads      = 10_000,
    mincomponents = 100    
    )
    df = load_rdata(; rdatafilename = datafilename)
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
    return standardized(df)
end

end # module OTUSampler
#/ End module
