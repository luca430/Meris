#= Simple script to investigate typical size of LEGO sets =#
#/ Packages
using DataFrames
using CSV

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "documentsize/lego/"
mkpath(DATADIR)

save = true
df = Meris.LegoSampler.computevocabsize(; returnsummary=true)
#~ Save
if save
    filename = "lego-documentsize.csv"
    CSV.write(DATADIR*filename, df)
end
