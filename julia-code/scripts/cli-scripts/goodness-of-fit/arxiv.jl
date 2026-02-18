#= Goodness of fit for arXiv papers =#
#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: arXivLoader, OhMyGoodness
using Meris: DATADIR

using JLD2

OUTDIR = DATADIR*"goodness-of-fit/arxiv/"
mkpath(OUTDIR)
FILENAME = "arxiv-candidatefits.jld2"

#/ Specify parameters
minsamples    = 30      #~ min. no. of articles per "class"
minreads      = 4_000   #~ min. no. of words per article
mincomponents = 10_000  #~ min. no. of distinct components per "class"
nε = 64

"Simple (private) function to load arXiv data"
function _load()
    #  note: data is filtered [using the parameters above]
    arxivdf = arXivLoader.load(
        stopwords=true, filterdata=true;
        minsamples=minsamples, minreads=minreads, mincomponents=mincomponents
    )
    #~ For each `domain`, fit CADs for each `sample_id` specifically, and store their parameters
    @transform!(arxivdf, :frequency = :counts ./ :nreads)
    return arxivdf
end

#/ Load and fit candidates [see `candidates.jl`]
arxivdf = _load()
fitdf = OhMyGoodness.fit_candidates(arxivdf; nε=nε)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf = fitdf)

# cols = [:GeneralizedPareto, :ParetoI, :ParetoIV, :TemperedPareto, :Gamma, :LogNormal]
# @rtransform! adf begin
#     :likely = begin
#         vals = Tuple(AsTable(cols))
#         cols[argmin(vals)]
#     end
# end

