#= Simple script to check Zipf's law in Wikitext-2 dataset =#
#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase
using JLD2

using Meris
#~ Specify directories to store data
DATADIR = Meris.DATADIR * "zipf/wikitext/"
mkpath(DATADIR)

save = true

DIRS = ["datasets/wikitext-2/", "datasets/wikitext-103/"]
FILENAMES = ["wikitext-2-raw.txt", "wikitext-103-raw.txt"]
sfilenames = ["wikitext-2", "wikitext-103"]

for i in eachindex(FILENAMES)
    #~ Load and make countmap
    wikitext = Meris.WikitextSampler.load_bagofwords(;
        FILENAME=FILENAMES[i],
        DIR=Meris.DATADIR*DIRS[i]
    )
    cm = countmap(wikitext)
    wikidf = DataFrame(component_id=collect(keys(cm)), count=collect(values(cm)))
    nreads = sum(wikidf[:,:count])
    @transform!(
        wikidf,
        :rank = tiedrank(:count, rev=true),
        :frequency = :count ./ nreads
    )

    if save
        filename = "$(sfilenames[i])-zipf.jld2"
        jldsave(DATADIR*filename; rank=wikidf.rank, frequency=wikidf.frequency)
    end
end
