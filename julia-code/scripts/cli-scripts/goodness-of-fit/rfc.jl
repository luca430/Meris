#= Candidate comparison and goodness-of-fit tests for RFC documents =#
#~ Parse command-line args
using Meris: MArgParse as Args
args = Args.parsegof()

#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: RFCLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

#/ Specify output directory and filename
OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "rfc-candidatefits.jld2"

@info "Fits and comparisons for RFC documents..."
#/ Load and fit candidates [see `candidates.jl`]
df = RFCLoader.load(; applyfilter=args["filter"], top=args["top"])
@transform!(df, :class = "rfc", :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    df, :class; testcandidate=:ParetoI, nε=args["numeps"], __computepvalue=false
)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf, aicdf=aicdf)
