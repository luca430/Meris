#= Simple script to verify expected scaling of new tokens with document size =#
#/ Packages
using DataFrames, DataFramesMeta
using JLD2

using Meris
DATADIR = Meris.DATADIR * "heaps/rfc/"
mkpath(DATADIR)

#~ Specify variables
save = true
nseeds = 144

#/ Compute the no. of distinct tokens (categories) in RFC documents
#~ Specify variables
save = true
mintokens = 512      #~ Min. no of tokens in text to be counted
maxfiles = 5000      #~ Max. no of parsed files
maxrows  = 1_000_000 #~ Max. no of rows allowed in DataFrame

#/ Compute the vocabularysize
heapsdf = Meris.RFCSampler.computevocabsize()
#~ Save
if save
    filename = "rfc-heaps.jld2"
    jldsave(DATADIR*filename; N=heapsdf.observationlength, V=heapsdf.vocabularysize)
end


