#= Simple script to compute the AFD for tokens in the RFC database =#
#/ Packages
using CSV, DataFrames

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "macro/afd/rfc/"
mkpath(DATADIR)

FILENAME = "z-values.csv"
TLFILENAME = "tl-stats.csv"

#~ Specify variables
save = true
mintokens = 1_000
maxfiles  = 2_000
maxrows   = 2^21

#/ Load data
rfcdf = Meris.RFCSampler.collect_rfcs(; mintokens=mintokens, maxfiles=maxfiles, maxrows=maxrows)
#/ Compute AFD and mean-variance
#~ note: here, :component_id is the relevant column
afddf = Meris.AFD.compute(rfcdf, :component_id; maxfrequency=1e-2)
tldf  = Meris.Taylor.compute(rfcdf, :component_id)

#/ Save
if save
    CSV.write(DATADIR*FILENAME, afddf)
    CSV.write(DATADIR*TLFILENAME, tldf)
end
