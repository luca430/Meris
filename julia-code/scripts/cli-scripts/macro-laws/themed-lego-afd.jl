#= Simple script to compute the AFD for LEGO pieces in a specific set =#
#/ Packages
using CSV, DataFrames

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "macro/afd/lego/"
mkpath(DATADIR)

FILENAME = "themed-z-values.csv"
TLFILENAME = "themed-tl-stats.csv"

#~ Specify variables
save = true
minquantity = 64          #~ Min. amount of LEGO pieces in a set
mindistinctpieces = 32    #~ Min. amount of distinct LEGO pieces in a set

#/ Load data
#~ note: if `theme_id=nothing`, selects the theme with the most sets [Star Wars]
theme_id = nothing
legodf, themedf = Meris.LegoSampler.parse_themes(;
    minquantity=minquantity, mindistinctpieces = mindistinctpieces,
    standardize=true, returnthemes=true
)
#~ subselect the theme with the most sets [Star Wars `sw`]
swlegodf = Meris.LegoSampler.select_theme(legodf, themedf)
#/ Compute AFD and mean-variance
#~ note: here, :component_id is the relevant column
afddf = Meris.AFD.compute(swlegodf, :component_id; minoccupancy=1e-1)
tldf  = Meris.Taylor.compute(swlegodf, :component_id, minoccupancy=1e-1)

#/ Save
if save
    CSV.write(DATADIR*FILENAME, afddf)
    CSV.write(DATADIR*TLFILENAME, tldf)
end
