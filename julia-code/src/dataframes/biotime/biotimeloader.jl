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
function load(
    ;
    DIR=BIOTIMEDIR * "raw-data/",
    FILENAME="biotime_v2_rawdata_2025",
    filterdata=true,  #~ Filter data
    minsamples=30,
    minreads=5_000,
    mincomponents=200,
    minsamplecomponents=70,
    verbose=false
)
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
        df = Meris.DataTools.df_filter(
            df,
            minreads=minreads,
            mincomponents=mincomponents,
            minsamplecomponents=minsamplecomponents,
            minsamples=minsamples
        )
    end

    return df
end

end # module BioTIMESampler
#/ End module

