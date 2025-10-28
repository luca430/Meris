#= Simple script to compute the AFD for LEGO pieces in the LEGO sets database =#
#/ Packages
using CSV, DataFrames

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "macro/afd/lego/"
mkpath(DATADIR)

FILENAME = "z-values.csv"
TLFILENAME = "tl-stats.csv"

#~ Specify variables
save = true
minquantity = 64          #~ Min. amount of LEGO pieces in a set
mindistinctpieces = 32    #~ Min. amount of distinct LEGO pieces in a set

#/ Load data
legodf = Meris.LegoSampler.filterlegos(;
    minquantity=minquantity, mindistinctpieces=mindistinctpieces
)
#/ Compute AFD and mean-variance
#~ note: here, :component_id is the relevant column
afddf = Meris.AFD.compute(legodf, :component_id; maxfrequency=1.)
tldf  = Meris.Taylor.compute(legodf, :component_id, maxfrequency=1.)

#/ Save
if save
    CSV.write(DATADIR*FILENAME, afddf)
    CSV.write(DATADIR*TLFILENAME, tldf)
end
