#= Simple script to check expected scaling of new categories with sample size in a Chinese book
   Somehow, as Chinese is more simple (character-wise), logarithmic scaling is expected
=#
#/ Packages
using JLD2

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "heap/books/chinese/"
mkpath(DATADIR)

#~ Specify variables
save = true
N = trunc.(Int, exp10.(range(2,log10.(700_000),16)))
nseeds = 144

#/ Compute the no. of categories V
#~ note: the function `count_categories` returns the mean no. of categories from a
#        Pitman-Yor process of size N
V = Meris.BookSampler.samplevocabsize(N; n=nseeds)
#~ Save
if save
    filename = "stone-vocabsize.jld2"
    jldsave(DATADIR*filename; V=V, N=N)
end
