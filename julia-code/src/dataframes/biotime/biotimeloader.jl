#= Module to load BioTIME dataset

=#
#/ Start module
module BioTIMELoader

#/ Packages
using ZipFile, CSV
using DataFrames, DataFramesMeta
using Meris

#/ Modules, directories
import ..DataTools: filterdata
import Meris.BIOTIMEDIR as BIOTIMEDIR

#################
### FUNCTIONS ###
import Meris.BIOTIMEDIR as BIOTIMEDIR

#################
### FUNCTIONS ###
function load(
    ;
    DIR = BIOTIMEDIR,
    FILENAME = "biotime_v2_rawdata_2025",
    fromzip       = true,
    fromparsed    = true,  #~ Load data that already has been parsed [by `savedata=true`]
    applyfilter   = true,  #~ Filter data
    savedata      = true,  #~ Store filtered data for easy retrieval
    minsamples    = 30,
    minreads      = 100_000,
    mincomponents = 100,    
    reorder       = true,
    top           = nothing,
    verbose       = false
    )
    #~ Load already parsed CSV if it exists
    FILTEREDFILENAME = "filtered_" * FILENAME * ".csv"
    if isfile(DIR*FILTEREDFILENAME) && fromparsed
        savedata = false    #~ No need to save again
        df = CSV.read(DIR*FILTEREDFILENAME, DataFrame)
    else
        if fromzip
            zip_path = DIR * "raw-data/$(FILENAME).zip"
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
            FILENAME = "raw-data/$(FILENAME).csv"
            df = CSV.read(DIR * FILENAME, DataFrame)
        end

        #~ Remove all "NA"
        subset!(df, All() .=> ByRow(!=("NA")))
        df.ABUNDANCE .= parse.(Float64, df.ABUNDANCE)
        df.sample_id .= string.(df.YEAR).* df.MONTH .* df.DAY
        @select!(df, :taxon, :STUDY_ID, :sample_id, :ID_SPECIES, :ABUNDANCE)
        df = @rename df begin
            :class = :STUDY_ID
            :component_id = :ID_SPECIES
            :counts = :ABUNDANCE
        end
        #~ Group and compute total no. of reads
        df = @chain df begin
            @groupby(:class, :sample_id)
            @combine(:component_id, :counts, :nreads=sum(:counts))
        end
        #~ Aggregate, as some samples were gathered on the same day
        df = @chain df begin
            @groupby(:class, :sample_id, :component_id)
            @combine(
                :counts = sum(:counts),
                :nreads = first(:nreads)
            )
        end
    end
        
    if applyfilter
        #~ filter data
        df = filterdata(
            df; minsamples=minsamples, minreads=minreads, mincomponents=mincomponents,
            reorder=reorder, top=top
        )
    end

    #~ Store for easy retrieval later
    (savedata) && (CSV.write(DIR * "$(FILTEREDFILENAME)", df))
    return df
end

end # module BioTIMESampler
#/ End module

