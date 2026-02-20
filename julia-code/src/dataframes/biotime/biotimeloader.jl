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
    DIR = BIOTIMEDIR * "raw-data/",
    FILENAME = "biotime_v2_rawdata_2025",
    fromzip       = true,
    fromparsed    = false, #~ Load data that already has been parsed [by `savedata=true`]
    filterdata    = true,  #~ Filter data
    savedata      = true,  #~ Store filtered data for easy retrieval
    minsamples    = 30,
    minreads      = 5_000,
    mincomponents = 100,
    verbose       = false
    )
    if !fromparsed
        if fromzip
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
            df = CSV.read(FILENAME, DataFrame)
            close(z)
        else
            #~ Load raw [extracted] CSV directly
            FILENAME = "$(FILENAME).csv"
            df = CSV.read(DIR * FILENAME, DataFrame)
        end

        #~ Remove all "NA"
        subset!(df, All() .=> ByRow(!=("NA")))
        df.ABUNDANCE .= parse.(Float64, df.ABUNDANCE)
        df.sample_id .= string.(df.YEAR).* df.MONTH .* df.DAY
        select!(df, :taxon, :STUDY_ID, :sample_id, :ID_SPECIES, :ABUNDANCE)
        rename!(df, :STUDY_ID => :class, :ID_SPECIES => :component_id, :ABUNDANCE => :counts)
        df = @chain df begin
            @groupby(:class, :sample_id)
            @combine(:component_id, :counts, :nreads=sum(:counts))
        end

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
    else
        #~ Load already parsed CSV
        savedata = false    #~ No need to save again
        FILENAME = "filtered_" * FILENAME * ".csv"
        df = CSV.read(DIR*FILENAME, DataFrame)
    end

    #~ Store for easy retrieval later
    (savedata) && (CSV.write(DIR * "filtered_$(FILENAME).csv", df))
    return df
end

end # module BioTIMESampler
#/ End module

