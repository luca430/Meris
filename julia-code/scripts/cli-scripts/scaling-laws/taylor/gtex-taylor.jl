#= Script to check Taylor's law in RFC documents =#
#/ Packages
using DataFrames, DataFramesMeta
using StatsBase
using JLD2

import Meris.DATADIR as DATADIR
OUTDIR = DATADIR * "taylor/gtex/"
mkpath(OUTDIR)

import Meris.GTExSampler as GTEx
import Meris.Taylor as Taylor

#/ Variables
tissues = ["BRAIN"]
# mincounts = 16
# mingenes = 128

#/ Load data
gtxdf = GTEx.load_gtex(; tissues=tissues)
nsamples = length(unique(gtxdf[!,:sample_id]))
tdf = Taylor.compute(gtxdf, :component_id; minoccupancy=32/nsamples)

#/ Save
filename = "gtex-taylor.jld2"
jldsave(
    OUTDIR*filename;
    mean=tdf.meanfrequency, var=tdf.varfrequency,     #~ Mean and variance
    omean=tdf.omeanfrequency, ovar=tdf.ovarfrequency, #~ Same, but taking into account occupancy
    occupancy=tdf.occupancy,                          #~ Raw occupancies
    varianceman=tdf.varm, variancevariance=tdf.vars,  #~ Sample variance of the mean and variance
    errorcorr=tdf.errorcorr, errorcov=tdf.errorcov    #~ Correlation and covariance of errors
)
