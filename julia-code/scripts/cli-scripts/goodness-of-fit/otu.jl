#= Candidate comparison and goodness-of-fit tests for OTU environments [EBI Metagenomics] =#
#~ Parse command-line args
using Meris: MArgParse as Args
args = Args.parsegof()

#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: OTULoader, OhMyGoodness
using Meris: DATADIR

using JLD2

#/ Specify output directory and filename
OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "otu-candidatefits.jld2"

@info "Fits and comparisons for EBI Metagenomics OTU count data..."
#/ Load and fit candidates
df = OTULoader.load(top=10)
@transform!(df, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    df, :class;
    testcandidate=:TemperedPareto, nε=args["numeps"]
)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf=fitdf, aicdf=aicdf)
