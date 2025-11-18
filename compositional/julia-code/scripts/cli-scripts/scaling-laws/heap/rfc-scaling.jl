#= Simple script to verify expected scaling of new tokens with document size =#
#/ Packages
using JLD2

using Meris
DATADIR = Meris.DATADIR * "heap/rfc/"
mkpath(DATADIR)

#~ Specify variables
save = true
nseeds = 128

#/ Compute the no. of distinct tokens (categories) in RFC documents
#@ TODO IMPLEMENT THIS
