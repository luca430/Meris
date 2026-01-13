#= Simple script to compute SAD for RFC documents =#
#/ Packages
using DataFrames, DataFramesMeta
using StatsBase
using JLD2

using Meris
#~ Specify directories to store data
DATADIR = Meris.DATADIR * "sad/rfc/"
mkpath(DATADIR)

#~ As there are /many/ RFC documents,
#  parse only a subset that fits in working memory
mintokens = 512       #~ Min. no of tokens in text to be counted
maxfiles = 128       #~ Max. no of parsed files
maxrows  = 1_000_000  #~ Max. no of rows allowed in DataFrame

#/ Load RFC data
rfcdf = Meris.RFCSampler.collect_rfcs(; mintokens=mintokens, maxfiles=maxfiles, maxrows=maxrows)
@transform!(rfcdf, :frequency = :counts ./ :nreads)


# fig = Figure(; size=(246, 246))
# ax = Axis(fig[1,1])

# observations = unique(rfcdf[:,:sample_id])
# for k in observations
#     logfreqs = log10.(rfcdf[rfcdf[!,:sample_id] .== k, :][!, :frequency])
#     minfreq, maxfreq = extrema(logfreqs)
#     binedges = collect(range(minfreq, maxfreq, 21))
#     fh = FHist.Hist1D(logfreqs; binedges=binedges, overflow=true)
#     fhx = FHist.bincenters(fh)[fh.bincounts .> 0.]
#     fhy = fh.bincounts[fh.bincounts .> 0.]
#     Makie.scatter!(ax, fhx, log10.(fhy), markersize=3)
# end
