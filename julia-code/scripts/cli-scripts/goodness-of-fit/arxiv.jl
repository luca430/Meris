#= Goodness of fit for arXiv papers =#
exit()
#~ Parse command-line args
using Meris: MArgParse as Args
args = Args.parsegof()

#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: arXivLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "arxiv-candidatefits.jld2"

#/ Specify parameters
minsamples    = 30      #~ min. no. of articles per "class"
minreads      = 4_000   #~ min. no. of words per article
mincomponents = 10_000  #~ min. no. of distinct components per "class"

@info "Fits and comparisons for arXiv articles..."
#/ Load and fit candidates [see `candidates.jl`]
df = arXivLoader.load(
    stopwords=true, applyfilter=args["filter"], top=args["top"],
    minsamples=minsamples, minreads=minreads, mincomponents=mincomponents
)
@transform!(df, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    df, :class; testcandidate=:ParetoI, nε=args["numeps"], __computepvalue=true
)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf = fitdf, aicdf = aicdf)

