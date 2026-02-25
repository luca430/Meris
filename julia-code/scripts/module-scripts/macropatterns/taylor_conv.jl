module TaylorConv

using Meris
using DataFrames, DataFramesMeta, StatsBase
using JLD2

function occupancy_scaling(df; occ_vec=collect(0.3:0.1:1.0), N=Int(1e5), save=true, filename="b_o_scaling.jld2")
    classes = unique(df.class)
    d = Dict()
    for class in classes
        println(class)
        sub = df[df.class .== class, :]
        b_vec = []
        for occ in occ_vec
            println("    ", occ)
            ds_counts = Meris.DataTools.downsample(sub; N=N, class=nothing)
            dsf_counts = Meris.DataTools.order_by_occ(ds_counts; occ=occ)
            x, y, fit = Meris.Taylor.compute2(dsf_counts; binned=true, nbins=30)
            push!(b_vec, fit.b)
        end
        d[class] = [occ_vec, b_vec]
    end
    
    (save) && (@save filename d)
    return d
end

function samplesize_scaling(df; size_vec=Int.(floor.(10 .^ collect(3:0.5:7.5))), occ=0.75, save=true, filename="b_N_scaling.jld2")
    classes = unique(df.class)
    d = Dict()
    for class in classes
        println(class)
        sub = df[df.class .== class, :]
        b_vec = []
        for N in size_vec
            println("    ", N)
            ds_counts = Meris.DataTools.downsample(sub; N=N, class=nothing)
            dsf_counts = Meris.DataTools.order_by_occ(ds_counts; occ=occ)
            x, y, fit = Meris.Taylor.compute2(dsf_counts; binned=true, nbins=30)
            push!(b_vec, fit.b)
        end
        d[class] = [size_vec, b_vec]
    end
    
    (save) && (@save filename d)
    return d
end

end # End module