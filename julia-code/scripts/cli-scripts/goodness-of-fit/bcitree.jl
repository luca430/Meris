#= Goodness of fit for Barro-Colorato Island tree counts =#
# exit()
#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: BCITreeLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

OUTDIR = DATADIR*"gof/"
mkpath(OUTDIR)
FILENAME = "bcitree-candidatefits.jld2"

#/ Specify parameters
# mincount = 50
nε = 100

#~ The BCI Tree data is a bit special, as there are quite few trees per quadrat
#  Therefore, fitting a heavy-tailed distribution, which requires typically some decades
#  of variation in (relative) counts is useless from the get-go. To this end, we simply
#  aggregate here and consider the entire island as a single `sample_id`.
#  @TODO LOAD THIS PROPERLY
treedf = BCITreeLoader.load()
treedf.class .= "eco-BCI"
@transform!(treedf, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    treedf, :class; testcandidate=:TemperedPareto, nε=100
)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf = fitdf, aicdf = aicdf)