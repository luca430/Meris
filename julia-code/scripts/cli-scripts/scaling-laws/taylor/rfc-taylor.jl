#= Script to check Taylor's law in RFC documents =#
#/ Packages
using DataFrames, DataFramesMeta
using StatsBase
using JLD2

import Meris.DATADIR as DATADIR
OUTDIR = DATADIR * "taylor/rfc/"
mkpath(OUTDIR)

import Meris.RFCSampler as RFC
import Meris.Taylor as Taylor

#/ Variables
mintokens = 512      #~ Min. no of tokens in text to be counted
maxfiles = 5000      #~ Max. no of parsed files
maxrows  = 1_000_000 #~ Max. no of rows allowed in DataFrame
#/ Load data
rfcdf = RFC.collect_rfcs(; mintokens=mintokens, maxfiles=maxfiles, maxrows=maxrows)
nsamples = length(unique(rfcdf[!,:sample_id]))
tdf = Taylor.compute(rfcdf, :component_id; minoccupancy=32/nsamples)

#/ Save
filename = "rfc-taylor.jld2"
jldsave(
    OUTDIR*filename;
    mean=tdf.meanfrequency, var=tdf.varfrequency,     #~ Mean and variance
    omean=tdf.omeanfrequency, ovar=tdf.ovarfrequency, #~ Same, but taking into account occupancy
    occupancy=tdf.occupancy,                          #~ Raw occupancies
    varianceman=tdf.varm, variancevariance=tdf.vars,  #~ Sample variance of the mean and variance
    errorcorr=tdf.errorcorr, errorcov=tdf.errorcov    #~ Correlation and covariance of errors
)
