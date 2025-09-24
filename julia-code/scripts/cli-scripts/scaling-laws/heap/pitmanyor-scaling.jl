#= Simple script to verify expected scaling of new categories with sample size =#
#/ Packages
using JLD2

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "heap/pitmanyor/"
mkpath(DATADIR)

#~ Specify variables
save = true
Nv = trunc.(Int, exp10.(range(2,5,16)))
αv = [0.0, 0.2, 0.5, 0.9]
θ  = 10.0

nseeds = 144

for α in αv
    #/ Compute the no. of categories V
    #~ note: the function `count_categories` returns the mean no. of categories from a
    #        Pitman-Yor process of size N
    V = Meris.PitmanYor.countvocabsize(Nv; n=nseeds, θ=θ, α=α)
    #~ Save
    if save
        αs = round(α, digits=1)
        filename = "pitmanyor-vocabsize_a$(αs).jld2"
        jldsave(DATADIR*filename; V=V, N=Nv, θ=θ)
    end
end
