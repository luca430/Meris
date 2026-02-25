#= Goodness of fit for Gutenberg book data =#
exit()
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

OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "gutenberg-candidatefits.jld2"

@info "Fits and comparisons for Gutenberg Project book data..."
#/ Load and fit candidates [see `candidates.jl`]
df = GutenbergLoader.load(; applyfilter=args["filter"], top=args["top"])
@transform!(df, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    df, :class; testcandidate=:ParetoI, nε=args["numeps"], __computepvalue=false
)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf, aicdf=aicdf)
