#= Goodness of fit for GTEx data =#
#~ Parse command-line args
using Meris: MArgParse as Args
args = Args.parsegof()
#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: GTExLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "gtex-candidatefits.jld2"

@info "Fits and comparisons for genetic (GTEx) data..."
#/ Load and fit candidates [see `candidates.jl`]
gtexdf = GTExLoader.load(top=10)
@transform!(gtexdf, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    gtexdf, :class; testcandidate=:ParetoI, nε=args["numeps"]
)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf, aicdf=aicdf)
