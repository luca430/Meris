#= Simple script to check Heaps' law in BCI tree data =#
#/ Packages
using DataFrames, DataFramesMeta
using StatsBase
using JLD2

using Meris
DATADIR = Meris.DATADIR * "heaps/bci.tree/"
mkpath(DATADIR)

#~ Specify variables
save = true

bcitreedf = Meris.BCITreeSampler.computevocabularysize()
#~ Save
if save
    filename = "bci.tree-heaps.jld2"
    jldsave(DATADIR*filename; N=bcitreedf.observationlength, V=bcitreedf.vocabularysize)
end
