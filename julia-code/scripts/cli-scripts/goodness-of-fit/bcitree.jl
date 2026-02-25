#= Goodness of fit for Barro-Colorato Island tree counts =#
#~ Parse command-line args
using Meris: MArgParse as Args
args = Args.parsegof()
#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: BCITreeLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "bcitrees-candidatefits.jld2"

@info "Fits and comparisons for Barro-Colorato Island data..."
#/ Load and fit candidates [see `candidates.jl`]
df = BCITreeLoader.load(top=10)
@transform!(df, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    df, :class;
    testcandidate=:ParetoIV, nε=args["numeps"]
)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf = fitdf, aicdf = aicdf)