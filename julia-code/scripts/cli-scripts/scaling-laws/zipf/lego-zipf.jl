#= Simple script to check Zipf's law in LEGO sets =#
#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase
using JLD2

using Meris
#~ Specify directories to store data
DATADIR = Meris.DATADIR * "zipf/lego/"
mkpath(DATADIR)
#~ Specify directories to find data
DIR = Meris.LEGODIR
FILENAME = "inventory_parts.csv"

#~ Specify variables
save = true
minquantity = 100      #~ Min. no of LEGO pieces in a set
mindistinctpieces = 30 #~ Min. no of *distinct* LEGO pieces in a set


legodf = Meris.LegoSampler.filterlegos(;
    minquantity=minquantity, mindistinctpieces=mindistinctpieces,
    DIR=DIR, FILENAME=FILENAME
)
rankdf = @chain legodf begin
    @by(:component_id, :totalcount = sum(:counts))
    @transform(:rank = tiedrank(:totalcount, rev=true))
    @transform(:frequency = :totalcount ./ sum(:totalcount))
end

#~ Save
if save
    filename = "lego-zipf.jld2"
    jldsave(DATADIR*filename; rank=rankdf.rank, frequency=rankdf.frequency)
end
