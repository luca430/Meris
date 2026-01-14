#= Simple script to verify expected scaling of new categories with sample size =#
#/ Packages
using DataFrames
using StatsBase
using JLD2

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "heaps/otu/"
mkpath(DATADIR)

#~ Specify variables
save = true
Nv = trunc.(Int, exp10.(range(1,log10.(12_000),25)))
nseeds = 144

environments = ["gut1", "gut2"]

for env in environments
    #/ Load data for that specific environment
    edf = Meris.OTUSampler.load_data(env)
    #~ filter data
    #!note: Defaults are given in `filter_data`
    fdf = Meris.OTUSampler.filter_data(edf)
    #/ Compute the no. of categories V in nseeds samples of sizes Nv=[N1,N2,...]
    #~ note: if `filterdf=true`, it filters the data, but we have done it before so that
    #        any values for filtering can be changed here instead of elsewhere, so set `false`
    V = Meris.OTUSampler.samplevocabsize(fdf, Nv; n=nseeds, filterdf=false)
    #~ Save
    if save
        filename = "otu-$(env)-heaps.jld2"
        jldsave(DATADIR*filename; N=Nv, V=dropdims(mean(V, dims=2), dims=2))
    end
end
