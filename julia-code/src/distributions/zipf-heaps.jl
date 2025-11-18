#= Simple module to compute Zipf's and Heap's law from standardized DataFrames =#
#/ Start module
module ZipfHeaps

using DataFrames, DataFramesMeta, StatsBase

#################
### FUNCTIONS ###
"Compute ranks and total counts from a standardized DataFrame. If 'class' is passed, the aggregation is done within the same class."
function zipf(df; class=nothing)
    if isnothing(class)
        sdf = @chain df begin
            @groupby(:component_id)
            @combine(:class_counts = sum(:counts))
            @transform(:ranks = tiedrank(-:class_counts))
        end
    else
        sdf = @chain df begin
            @groupby(class, :component_id)
            @combine(:class_counts = sum(:counts))
            @groupby(class)
            @combine(
                :component_id,
                :ranks = tiedrank(-:class_counts),
                :class_counts
            )
        end
    end
    return sdf
end


end # module ParetoLike
#/ End module
