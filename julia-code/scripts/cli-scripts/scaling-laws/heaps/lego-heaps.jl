#= Simple script to verify expected scaling of new categories with sample size =#
#/ Packages
using StatsBase
using JLD2

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "heaps/lego/"
mkpath(DATADIR)

#~ Specify variables
save = true
Nv = trunc.(Int, exp10.(range(2,5,16)))
nseeds = 256

#/ Compute the no. of categories V in nseeds samples of sizes Nv=[N1,N2,...]
V = Meris.LegoSampler.samplevocabsize(Nv; n=nseeds)
#~ Save
if save
    filename = "lego-heaps.jld2"
    jldsave(DATADIR*filename; V=dropdims(mean(V, dims=2), dims=2), N=Nv)
end
