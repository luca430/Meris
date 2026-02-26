module Taylor

using Meris
using DataFrames, DataFramesMeta, StatsBase
using JLD2

function compute(
        df;
        save = false,
        filename = "tl.jld2"
    )

    #/ Get specific taylor df. Do this one class at a time to preserve memory in the computation
    classes = unique(df.class)
    dfs = []
    for class in classes
        println(class)
        sub = df[df.class .== class, :]
        _tldf = Meris.Taylor.compute(sub, :component_id)
        _tldf.class .= class
        push!(dfs, _tldf) 
    end

    tldf = vcat(dfs...)

    (save) && (@save filename tldf)

    return tldf
end

end # End module