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
function load(; DIR=BIOTIMEDIR * "raw-data/", fromzip=true, verbose=false)
    if fromzip
        zip_path = DIR * "biotime_v2_rawdata_2025.zip"
        z = ZipFile.Reader(zip_path)

        # List files inside zip (if verbose=true)
        if verbose
            for f in z.files
                println(f.name)
            end
        end

        # Find the main CSV file (example name, adjust if needed)
        file = only(filter(f -> endswith(f.name, ".csv"), z.files))
        df = CSV.read(file, DataFrame)
        close(z)
    else
        #~ Load CSV directly
        file = DIR * "biotime_v2_rawdata_2025.csv"
        df = CSV.read(file, DataFrame)
    end

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

    return df
end

end # module BioTIMESampler
#/ End module

