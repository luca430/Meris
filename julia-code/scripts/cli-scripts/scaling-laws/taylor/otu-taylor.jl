#= Script to check Taylor's law for OTUs =#
#/ Packages
using DataFrames, DataFramesMeta
using StatsBase
using JLD2

import Meris.DATADIR as DATADIR
OUTDIR = DATADIR * "taylor/otu/"
mkpath(OUTDIR)

import Meris.OTUSampler as OTU
import Meris.Taylor as Taylor

#/ Variables
env = "gut1"
minsamples = 30
minreads = 10_000
mincounts = 16
minspecies = 24

#/ Load data
otudf = OTU.filter_data(
    OTU.load_data(env);
    minsamples=minsamples, minreads=minreads, mincounts=mincounts, minspecies=minspecies
)
nsamples = length(unique(otudf[!,:sample_id]))
tdf = Taylor.compute(otudf, :otu_id; minoccupancy=0.05)

#/ Save
filename = "otu-taylor.jld2"
jldsave(
    OUTDIR*filename;
    mean=tdf.meanfrequency, var=tdf.varfrequency,     #~ Mean and variance
    omean=tdf.omeanfrequency, ovar=tdf.ovarfrequency, #~ Same, but taking into account occupancy
    occupancy=tdf.occupancy,                          #~ Raw occupancies
    varianceman=tdf.varm, variancevariance=tdf.vars,  #~ Sample variance of the mean and variance
    errorcorr=tdf.errorcorr, errorcov=tdf.errorcov    #~ Correlation and covariance of errors
)
