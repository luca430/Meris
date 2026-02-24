#= Goodness of fit for LEGO sets =#
#~ Parse command-line args
using Meris: MArgParse as Args
args = Args.parsegof()

#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: LegoLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

#/ Specify output directory and filename
OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "lego-candidatefits.jld2"

@info "Fits and comparisons for LEGO set data..."
#/ Load and fit candidates [see `candidates.jl`]
df = LegoLoader.load(; filterdata=true, aggregate=true)
@transform!(df, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    df,:class; testcandidate=:ParetoIV, nε=args["numeps"]
)
exit()
