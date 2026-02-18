#= Goodness of fit for Barro-Colorato Island tree counts =#
#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: GTExLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

OUTDIR = DATADIR*"goodness-of-fit/gtex/"
mkpath(OUTDIR)
FILENAME = "gtex-candidatefits.jld2"

#/ Specify parameters
nε = 64

gtexdf = GTExLoader.load()
@transform!(gtexdf, :frequency = :counts ./ :nreads)
fitdf = OhMyGoodness.fit_candidates(gtexdf, :class; nε=nε)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf)
