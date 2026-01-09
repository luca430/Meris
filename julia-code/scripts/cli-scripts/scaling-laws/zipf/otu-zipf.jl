#= Simple script to check Zipf's law in LEGO sets =#
#/ Packages
using CSV, DataFrames, DataFramesMeta
using StatsBase
using JLD2

using Meris
#~ Specify directories to store data
DATADIR = Meris.DATADIR * "zipf/otu/"
mkpath(DATADIR)
#~ Specify directories to find data
DIR = Meris.OTUDIR

#~ Specify variables
save = true

environments = ["gut1", "gut2", "seawater"]
for env in environments
    #/ Load data for that specific environment
    edf = Meris.OTUSampler.load_data(env)
    #~ filter data
    #!note: Defaults are given in `filter_data`
    fdf = Meris.OTUSampler.filter_data(edf)

    rankdf = @chain fdf begin
        @by(:otu_id, :totalcount = sum(:counts))
        @transform(:rank = tiedrank(:totalcount, rev=true))
        @transform(:frequency = :totalcount ./ sum(:totalcount))
    end

    #~ Save
    if save
        filename = "otu-$(env)-zipf.jld2"
        jldsave(DATADIR*filename; rank=rankdf.rank, frequency=rankdf.frequency)
    end
end
