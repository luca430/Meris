#= Goodness of fit for Gutenberg book data =#
@info "Fits and comparisons for Gutenberg Project book data..."
#~ Parse command-line args
using Meris: MArgParse as Args
args = Args.parsegof()
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

#/ Load and fit candidates [see `candidates.jl`]
df = GutenbergLoader.load()
@transform!(df, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(df, :class; nε=args["numeps"])

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf, aicdf=aicdf)
