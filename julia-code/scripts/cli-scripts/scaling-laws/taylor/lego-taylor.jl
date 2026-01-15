#= Script to check Taylor's law in Lego sets =#
#/ Packages
using DataFrames, DataFramesMeta
using StatsBase
using JLD2

import Meris.DATADIR as DATADIR
OUTDIR = DATADIR * "taylor/lego/"
mkpath(OUTDIR)

import Meris.LegoSampler as Lego
import Meris.Taylor as Taylor

#/ Variables
minquantity = 100
mindistinctpieces = 50
#/ Load data
legodf = Lego.filterlegos(; minquantity=minquantity, mindistinctpieces=mindistinctpieces)
nsamples = length(unique(legodf[!,:sample_id]))
tdf = Taylor.compute(legodf, :component_id; minoccupancy=0.)

#/ Save
filename = "lego-taylor.jld2"
jldsave(
    OUTDIR*filename;
    mean=tdf.meanfrequency, var=tdf.varfrequency,     #~ Mean and variance
    omean=tdf.omeanfrequency, ovar=tdf.ovarfrequency, #~ Same, but taking into account occupancy
    occupancy=tdf.occupancy,                          #~ Raw occupancies
    varianceman=tdf.varm, variancevariance=tdf.vars,  #~ Sample variance of the mean and variance
    errorcorr=tdf.errorcorr, errorcov=tdf.errorcov    #~ Correlation and covariance of errors
)
