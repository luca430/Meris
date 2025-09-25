#= Simple script to verify expected scaling of new categories with sample size =#
#/ Packages
using JLD2

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "heap/lego/"
mkpath(DATADIR)

#~ Specify variables
save = true
Nv = trunc.(Int, exp10.(range(2,5,16)))
nseeds = 144 ÷ 2

#/ Compute the no. of categories V in nseeds samples of sizes Nv=[N1,N2,...]
V = Meris.LegoSampler.samplevocabsize(Nv; n=nseeds)
#~ Save
if save
    filename = "lego-vocabsize.jld2"
    jldsave(DATADIR*filename; V=V, N=Nv)
end
