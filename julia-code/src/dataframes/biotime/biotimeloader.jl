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
function load(
    ;
    DIR=BIOTIMEDIR * "raw-data/",
    FILENAME="biotime_v2_rawdata_2025",
    filterdata=true,  #~ Filter data
    minsamples=30,
    minreads=5_000,
    mincomponents=200,
    top=50,
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
    df = CSV.read(FILENAME, DataFrame, select=[:ABUNDANCE, :SAMPLE_DESC, :taxon, :STUDY_ID, :ID_SPECIES])
    close(z)
    
    #~ Remove all "NA" and rename columns to match our standard format
    subset!(df, All() .=> ByRow(!=("NA")))
    abundance_parsed = tryparse.(Int64, string.(df.ABUNDANCE))
    valid_abundance = .!isnothing.(abundance_parsed)
    if verbose && any(.!valid_abundance)
        println("Dropped $(count(.!valid_abundance)) rows with non-integer ABUNDANCE values")
    end
    df = df[valid_abundance, :]
    df.ABUNDANCE = Int64[x::Int64 for x in abundance_parsed[valid_abundance]]
    select!(df, :taxon, :STUDY_ID, :SAMPLE_DESC, :ID_SPECIES, :ABUNDANCE)
    rename!(df, :STUDY_ID => :class, :SAMPLE_DESC => :sample_id, :ID_SPECIES => :component_id, :ABUNDANCE => :counts)

    #~ Groupby to get nreads per sample
    df = @chain df begin
        @groupby(:class, :sample_id)
        @combine(:component_id, :counts, :nreads=sum(:counts))
    end
    
    if filterdata
        df = Meris.DataTools.filterdata(
            df,
            minreads=minreads,
            mincomponents=mincomponents,
            minsamples=minsamples,
            top=top
        )
    end

    return df
end

end # module BioTIMESampler
#/ End module
