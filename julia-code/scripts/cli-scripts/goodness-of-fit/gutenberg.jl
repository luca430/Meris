#= Goodness of fit for Barro-Colorato Island tree counts =#
#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: GutenbergLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

OUTDIR = DATADIR*"goodness-of-fit/gutenberg/"
mkpath(OUTDIR)
FILENAME = "gutenberg-candidatefits.jld2"

#/ Specify parameters
nε = 64

bookdf = GutenbergLoader.load()
@transform!(bookdf, :frequency = :counts ./ :nreads)
fitdf = OhMyGoodness.fit_candidates(bookdf, :class; nε=nε)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf)
