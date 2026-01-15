#= Module to investigate macro(-ecological) laws in genes data from the GTEx project

=#
#/ Start module
module GTExSampler

#/ Packages
using DataFrames, DataFramesMeta
using Random, StatsBase
using CodecZlib, CSV

#/ Modules, directories
import Meris.GTEXDIR as GTEXDIR

const TISSUES_MAP = Dict(
                        "BRAIN" => "brain_cortex",
                        "HEART" => "heart_left_ventricle",
                        "LIVER" => "liver",
                        "PANCREAS" => "pancreas"
                    )

#################
### FUNCTIONS ###
function load_gtex(; tissues=["BRAIN"])
    # Load file with metadata info to get nreads for each sample
    fname = GTEXDIR * "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
    df_reads = CSV.read(fname, DataFrame; delim='\t')
    df_reads = select(df_reads, :SAMPID, :SMRDTTL)
    rename!(df_reads, :SAMPID => :sample_id, :SMRDTTL => :nreads)
    
    # Loop over different datasets to make a unique DataFrame
    all_df = []
    for tissue in tissues
        fname = GTEXDIR * "gene_reads_v10_$(TISSUES_MAP[tissue]).gct.gz"
        
        # Open a decompression stream over the gz file
        gz = GzipDecompressorStream(open(fname, "r"))
        
        # Read into a DataFrame
        df_raw = CSV.File(
            gz;
            delim='\t',
            header=3,
            strict=false,    # allow rows with varying column counts
        ) |> DataFrame
        
        close(gz)
        
        select!(df_raw, Not([:Description]))
        df = stack(df_raw, Not(:Name); variable_name = :sample_id, value_name = :counts)
        filter!(:counts => x -> x .> 0, df)
        rename!(df, :Name => :component_id)
        df.class .= tissue
        select!(df, :class, :component_id, :sample_id, :counts)
        df = leftjoin(df, df_reads, on = :sample_id)
        push!(all_df, df)
    end
    
    return vcat(all_df...)
end

end # module arXivSampler
#/ End module

