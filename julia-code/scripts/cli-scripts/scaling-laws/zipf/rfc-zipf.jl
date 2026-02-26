#= Simple script to check Zipf's law in LEGO sets =#
#/ Packages
using DataFrames, DataFramesMeta
using StatsBase
using JLD2

using Meris
DATADIR = Meris.DATADIR * "zipf/rfc/"
mkpath(DATADIR)

#~ Specify variables
save = true
mintokens = 512      #~ Min. no of tokens in text to be counted
maxfiles = 5000      #~ Max. no of parsed files
maxrows  = 1_000_000 #~ Max. no of rows allowed in DataFrame

rfcdf = Meris.RFCSampler.collect_rfcs()
#~ Construct bag of trees
bagofwords = @by(rfcdf, :component_id, :totalcount = sum(:counts))
rankdf = @transform(bagofwords,
    :rank=tiedrank(:totalcount, rev=true),
    :frequency = :totalcount ./ sum(:totalcount)
)

#~ Save
if save
    filename = "rfc-zipf.jld2"
    jldsave(DATADIR*filename; rank=rankdf.rank, frequency=rankdf.frequency)
end
