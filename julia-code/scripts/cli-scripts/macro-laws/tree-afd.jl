#= Simple script to compute the AFD for tokens in the RFC database =#
#/ Packages
using CSV, DataFrames

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "macro/afd/bci.tree/"
mkpath(DATADIR)

FILENAME = "z-values.csv"
TLFILENAME = "tl-stats.csv"

#~ Specify variables
save = true
mincount = 1

#/ Load data
treedf = Meris.BCITreeSampler.load_treedata(;
    filename="bci.tree8.rdata", mincount=mincount, joinquadrats=true
)
#/ Compute AFD and mean-variance
#~ note: here, :component_id is the relevant column
afddf = Meris.AFD.compute(treedf, :component_id)
tldf  = Meris.Taylor.compute(treedf, :component_id)

#/ Save
if save
    CSV.write(DATADIR*FILENAME, afddf)
    CSV.write(DATADIR*TLFILENAME, tldf)
end
