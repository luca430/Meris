#= Simple script to check Zipf's law in LEGO sets =#
#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase
using JLD2

using Meris
DATADIR = Meris.DATADIR * "zipf/bci.tree/"
mkpath(DATADIR)

#~ Specify variables
save = true

bcitreedf = Meris.BCITreeSampler.load_treedata(; joinquadrats=false)
#~ Construct bag of trees
bagoftrees = @by(bcitreedf, :component_id, :totalcount = sum(:counts))
rankdf = @transform(bagoftrees,
    :rank=tiedrank(:totalcount, rev=true),
    :frequency = :totalcount ./ sum(:totalcount)
)

#~ Save
if save
    filename = "bcitree-zipf.jld2"
    jldsave(DATADIR*filename; rank=rankdf.rank, frequency=rankdf.frequency)
end
