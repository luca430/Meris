#= Goodness of fit for BioTIME species counts =#
#~ Parse command-line args
using Meris: MArgParse as Args
args = Args.parsegof()
#~ Load some packages
using DataFrames, DataFramesMeta
using Distributions
#~ Specify some Meris modules and directories
using Meris: BioTIMELoader, OhMyGoodness
using Meris: DATADIR

using JLD2

OUTDIR = DATADIR*"goodness-of-fit/"
mkpath(OUTDIR)
FILENAME = "biotime-candidatefits.jld2"

@info "Fits and comparisons for BioTIME data..."
#/ Load and fit candidates [see `candidates.jl`]
df = BioTIMELoader.load(top=50)
#/ BioTIME has a lot of samples with very few species so we manually drop samples with less than 100 species
df = @chain df begin
    @groupby(:class, :sample_id)
    @combine(:component_id, :counts, :nreads, :ncomponents=length(unique(:component_id)))
end
df = df[df.ncomponents .> 100, :]
@transform!(df, :frequency = :counts ./ :nreads)
fitdf, aicdf = OhMyGoodness.fit_candidates(
    df, :class;
    testcandidate=:ParetoIV, nε=args["numeps"]
)

#/ Store
jldsave(OUTDIR*FILENAME; fitdf = fitdf, aicdf = aicdf)

