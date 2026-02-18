module Meris

## GLOBAL CONSTANTS
const minsamples::Int = 30
const mincomponents::Int = 100

## DIRECTORIES
const DATADIR = normpath(joinpath(@__DIR__, "..", "data/"))
const FIGDIR = normpath(joinpath(@__DIR__, "..", "figures/"))

const ARXIVDIR = DATADIR * "datasets/arxiv/"
const TREEDIR = DATADIR * "datasets/bci.tree/"
const BIOTIMEDIR = DATADIR * "datasets/biotime/"
const BRIGHTKITEDIR = DATADIR * "datasets/brightkite/"
const EMAILDIR = DATADIR * "datasets/emails/"
const FINANCEDIR = DATADIR * "datasets/finance/"
const GAIADIR = DATADIR * "datasets/gaia/"
const GOWALLADIR = DATADIR * "datasets/gowalla/"
const GTEXDIR = DATADIR * "datasets/gtex/"
const GBIFDIR = DATADIR * "datasets/gbif/"
const GUTENBERGDIR = DATADIR * "datasets/gutenberg/"
const LEGODIR = DATADIR * "datasets/lego/"
const OTUDIR = DATADIR * "datasets/otu/"
const RFCDIR = DATADIR * "datasets/rfc/"
const TARADIR = DATADIR * "datasets/taraocean/"

## SUBMODULES
#~ Data handlers, loaders, and samplers
include("dataframes/arxiv/arxivloader.jl")
include("dataframes/bci.tree/bcitreeloader.jl")
include("dataframes/biotime/biotimeloader.jl")
include("dataframes/brightkite/brightkiteloader.jl")
include("dataframes/emails/emailsloader.jl")
include("dataframes/finance/financeloader.jl")
include("dataframes/gaia/gaialoader.jl")
include("dataframes/gowalla/gowallaloader.jl")
include("dataframes/GTEx/gtexloader.jl")
include("dataframes/gutenberg/gutenbergloader.jl")
include("dataframes/lego/legoloader.jl")
include("dataframes/otu/otuloader.jl")
include("dataframes/rfc/rfcloader.jl")

#~ Processes
include("processes/dirichlet.jl")
include("processes/pitman-yor.jl")

# include("fits/power-law.jl")
# include("fits/double-power-law.jl")
# include("fits/mle.jl")
# include("fits/heapsmodel.jl")

#~ Utilities
include("dataframes/afd.jl")
include("dataframes/taylor.jl")
include("dataframes/datatools.jl")
include("distributions/distributions.jl")
include("distributions/lr_distributions.jl")

#~ Fitting
include("fits/candidate-distributions.jl")
include("fits/goodness-of-fit.jl")
include("fits/straight-line.jl")

end # module Meris
