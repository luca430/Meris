#= Simple script to verify expected scaling of new categories with sample size =#
#/ Packages
using StatsBase
using JLD2

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "heaps/pitmanyor/"
mkpath(DATADIR)

#~ Specify variables
save = true
Nv = trunc.(Int, exp10.(range(1,4,11)))
# αv = [0.1, 0.5, 0.9]
αv = [0.5, 0.9]
θ  = 10.0

nseeds = 32

for α in αv
    #/ Compute the no. of categories V
    #~ note: the function `count_categories` returns the mean no. of categories from a
    #        Pitman-Yor process of size N
    V = Meris.PitmanYor.countvocabsize(Nv; n=nseeds, θ=θ, α=α)
    #~ Save
    if save
        αs = round(α, digits=1)
        filename = "pitmanyor-heaps_a$(αs).jld2"
        jldsave(DATADIR*filename; N=Nv, V=dropdims(mean(V, dims=2), dims=2), θ=θ)
    end
end
