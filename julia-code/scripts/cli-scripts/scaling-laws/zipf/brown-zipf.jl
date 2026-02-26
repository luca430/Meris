#= Simple script to check Zipf's law for the Brown corpus =#
using DataFrames, DataFramesMeta
using StatsBase
using JLD2

using Meris
#~ Specify directories to store data
DATADIR = Meris.DATADIR * "zipf/brown/"
mkpath(DATADIR)

save = true

bagofwords = Meris.WordSampler.load_bagofwords()
cm = countmap(bagofwords)
df = DataFrame(component_id=collect(keys(cm)), count=collect(values(cm)))

nreads = sum(df[:,:count])
@transform!(
    df,
    :rank = tiedrank(:count, rev=true),
    :frequency = :count ./ nreads
)

if save
    filename = "brown-zipf.jld2"
    jldsave(DATADIR*filename; rank=df.rank, frequency=df.frequency)
end
