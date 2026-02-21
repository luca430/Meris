#= Module to load BioTIME dataset

=#
#/ Start module
module BioTIMELoader

#/ Packages
using ZipFile, CSV
using DataFrames, DataFramesMeta
using Meris

#/ Modules, directories
import Meris.BIOTIMEDIR as BIOTIMEDIR

#################
### FUNCTIONS ###
import Meris.BIOTIMEDIR as BIOTIMEDIR

#################
### FUNCTIONS ###
function load(
    ;
    DIR=BIOTIMEDIR * "raw-data/",
    FILENAME="biotime_v2_rawdata_2025",
    fromparsed=false, #~ Load data that already has been parsed [by `savedata=true`]
    filterdata=true,  #~ Filter data
    savedata=true,  #~ Store filtered data for easy retrieval
    minsamples=30,
    minreads=5_000,
    mincomponents=200,
    minsamplecomponents=70,
    verbose=false
)
    if !fromparsed
        zip_path = DIR * "$(FILENAME).zip"
        z = ZipFile.Reader(zip_path)

        # List files inside zip (if verbose=true)
        if verbose
            for f in z.files
                println(f.name)
            end
        end

        # Find the main CSV file
        FILENAME = only(filter(f -> endswith(f.name, ".csv"), z.files))
        df = CSV.read(FILENAME, DataFrame, select=[:ABUNDANCE, :YEAR, :MONTH, :DAY, :taxon, :STUDY_ID, :ID_SPECIES])
        close(z)

        #~ Remove all "NA"
        subset!(df, All() .=> ByRow(!=("NA")))
        df.ABUNDANCE .= parse.(Float64, df.ABUNDANCE)
        df.sample_id .= string.(df.YEAR) .* df.MONTH .* df.DAY
        select!(df, :taxon, :STUDY_ID, :sample_id, :ID_SPECIES, :ABUNDANCE)
        rename!(df, :STUDY_ID => :class, :ID_SPECIES => :component_id, :ABUNDANCE => :counts)
        df = @chain df begin
            @groupby(:class, :sample_id)
            @combine(:component_id, :counts, :nreads = sum(:counts))
        end

        if filterdata
            #~ filter data
            @subset!(df, :nreads .> minreads)
            df = @chain df begin
                @groupby(:class)
                @combine(:sample_id, :component_id, :counts, :nreads, :ncomponents = length(unique(:component_id)))
                @subset(:ncomponents .> mincomponents)
                @groupby(:class, :sample_id)
                @combine(:sample_id, :component_id, :counts, :nreads, :ncomponentspersample = length(unique(:component_id)))
                @subset(:ncomponentspersample .> minsamplecomponents)
                @groupby(:class)
                @combine(:sample_id, :component_id, :counts, :nreads, :nsamples = length(unique(:sample_id)))
                @subset(:nsamples .> minsamples)
            end
        end
    else
        #~ Load already parsed CSV
        savedata = false    #~ No need to save again
        FILENAME = "filtered_" * FILENAME * ".csv"
        df = CSV.read(DIR * FILENAME, DataFrame)
    end

    #~ Store for easy retrieval later
    (savedata) && (CSV.write(DIR * "$(FILTEREDFILENAME)", df))
    return df
end

end # module BioTIMESampler
#/ End module

