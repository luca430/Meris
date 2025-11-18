#= Simple script to investigate typical size of arXiv papers =#
#/ Packages
using DataFrames
using DataFramesMeta
using CSV

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "documentsize/arxiv/"
mkpath(DATADIR)

#~ Specify variables
CATEGORIES = readlines(Meris.ARXIVDIR*"categories.txt")
save = true

for CATEGORY in CATEGORIES
    arxivdf = Meris.arXivSampler.load_papers(; CATEGORY=CATEGORY)
    sdf = Meris.arXivSampler.summarize(arxivdf)

    #~ Save
    if save
        filename = "arxiv-$(CATEGORY)-documentsize.csv"
        CSV.write(DATADIR*filename, sdf)
    end    
end
