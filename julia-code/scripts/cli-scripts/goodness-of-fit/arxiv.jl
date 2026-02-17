#= Goodness of fit for arXiv papers =#
using DataFrames, DataFramesMeta

using Meris: arXivLoader, MDistributions

#/ Specify parameters
minsamples    = 30      #~ min. no. of articles per "class"
minreads      = 4_000   #~ min. no. of words per article
mincomponents = 10_000  #~ min. no. of distinct components per "class"

#/ Load
#  note: data is filtered [using the parameters above]
arxivdf = Meris.arXivLoader.load(
    stopwords=true, filterdata=true;
    minsamples=minsamples, minreads=minreads, mincomponents=mincomponents
)

#~ For each domain, fit CADs for each `sample_id` specifically, and store them
