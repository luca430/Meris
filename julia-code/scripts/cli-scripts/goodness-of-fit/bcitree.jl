#= Goodness of fit for Barro-Colorato Island tree counts =#
#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: BCITreeLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

OUTDIR = DATADIR*"goodness-of-fit/bcitree/"
mkpath(OUTDIR)
FILENAME = "bcitree-candidatefits.jld2"

#/ Specify parameters
mincount = 50
nε = 64

#~ The BCI Tree data is a bit special, as there are quite few trees per quadrat
#  Therefore, fitting a heavy-tailed distribution, which requires typically some decades
#  of variation in (relative) counts is useless from the get-go. To this end, we simply
#  aggregate here and consider the entire island as a single `sample_id`.
#  @TODO LOAD THIS PROPERLY
treedf = BCITreeLoader.load(; mincount=mincount, joinquadrats=false)
@transform!(treedf, :island = "bci", :sample_id = 1, :frequency = :counts ./ :nreads)
fitdf = OhMyGoodness.fit_candidates(treedf, :island; nε=nε)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf)
