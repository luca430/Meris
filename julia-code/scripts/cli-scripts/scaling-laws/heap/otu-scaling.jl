#= Simple script to verify expected scaling of new categories with sample size =#
#/ Packages
using CSV, DataFrames
using JLD2

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "heap/otu/"
mkpath(DATADIR)

#~ Specify variables
save = true
Nv = trunc.(Int, exp10.(range(2,log10.(12_000),14)))
nseeds = 144

environments = ["gut1"]

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
        filename = "otu-$(env)-vocabsize.jld2"
        jldsave(DATADIR*filename; V=V, N=Nv)
    end
end
