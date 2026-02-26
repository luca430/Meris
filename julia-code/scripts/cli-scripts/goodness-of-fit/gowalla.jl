#= Candidate comparison and goodness-of-fit tests for the Gowalla project =#
exit()
#~ Parse command-line args
# using Meris: MArgParse as Args
# args = Args.parsegof()

#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: GowallaLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

#/ Specify output directory and filename
OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "gowalla-candidatefits.jld2"

@info "Fits and comparisons for the Gowalla project..."
#/ Load and fit candidates [see `candidates.jl`]
df = GowallaLoader.load(; applyfilter=args["filter"], top=args["top"])
@transform!(df, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    df, :class; testcandidate=:ParetoI, nε=args["numeps"], __computepvalue=true
)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf, aicdf)
