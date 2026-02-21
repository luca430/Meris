#= Candidate comparison and goodness-of-fit tests for the stock volume data =#
#~ Parse command-line args
using Meris: MArgParse as Args
args = Args.parsegof()

#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: FinanceLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

#/ Specify output directory and filename
OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "finance-candidatefits.jld2"

@info "Fits and comparisons for Yahoo Finance stock volume data..."
#/ Load and fit candidates [see `candidates.jl`]
df = FinanceLoader.load(; filterdata=true)
@transform!(df, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(df,:class; testcandidate=:ParetoI, nε=args["numeps"])

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf, aicdf=aicdf)
