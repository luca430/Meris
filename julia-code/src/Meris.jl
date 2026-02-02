module Meris

## DIRECTORIES
const DATADIR   = normpath(joinpath(@__DIR__, "..", "data/"))
const FIGDIR = normpath(joinpath(@__DIR__, "..", "figures/"))

const ARXIVDIR = DATADIR * "datasets/arxiv/"
const BOOKDIR   = DATADIR * "datasets/books/"
const GTEXDIR   = DATADIR * "datasets/gtex/"
const GUTENBERGDIR = DATADIR * "datasets/gutenberg/"
const LEGODIR   = DATADIR * "datasets/lego/"
const OTUDIR    = DATADIR * "datasets/otu/"
const CORPUSDIR = DATADIR * "datasets/corpus/"
const RFCDIR    = DATADIR * "datasets/rfc/"
const TREEDIR   = DATADIR * "datasets/bci.tree/"
const TARADIR = DATADIR * "datasets/taraocean/"
const PATENTDIR = DATADIR * "datasets/patents/"
const WIKIDIR   = DATADIR * "datasets/wikitext/"
const WIKI2DIR  = DATADIR * "datasets/wikitext-2/"
const WIKI103DIR= DATADIR * "datasets/wikitext-103/"
const YAHOODIR = DATADIR * "datasets/finance/"

## SUBMODULES
#~ Data handlers, loaders, and samplers 
# include("dataframes/books/booksampler.jl")
include("dataframes/gtex/gtexsampler.jl")
# include("dataframes/books/wordsampler.jl")
include("dataframes/wikitext/wikitextsampler.jl")
include("dataframes/lego/legosampler.jl")
# include("dataframes/otu/otusampler.jl")
# include("dataframes/arxiv/arxivsampler.jl")
include("dataframes/rfc/rfcsampler.jl")

## SUBMODULES
#~ Data handlers, loaders, and samplers
include("dataframes/arxiv/arxivloader.jl")
include("dataframes/bci.tree/bcitreesampler.jl")
include("dataframes/gaia/gaialoader.jl")
include("dataframes/gowalla/gowallaloader.jl")
include("dataframes/GTEx/gtexloader.jl")
# include("dataframes/gutenberg/gutenbergloader.jl")
include("dataframes/lego/legosampler.jl")
include("dataframes/otu/otuloader.jl")
include("dataframes/rfc/rfcsampler.jl")
# include("dataframes/tara/taraloader.jl")

#~ Processes
include("processes/dirichlet.jl")
include("processes/pitman-yor.jl")

#~ Fitting
include("fits/goodness-of-fit.jl")
include("fits/straight-line.jl")

#~ Utilities
include("dataframes/afd.jl")
include("dataframes/taylor.jl")
include("dataframes/datatools.jl")

include("distributions/distributions.jl")

end # module Meris
