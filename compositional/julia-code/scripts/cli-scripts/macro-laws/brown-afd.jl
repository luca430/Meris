#= Simple script to compute the AFD for the Brown corpus =#
#/ Packages
using CSV, DataFrames
using Distributions
using Random
using JLD2

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "macro/afd/"
mkpath(DATADIR)

FILENAME = "zvalues.csv"
TLFILENAME = "stats.csv"

#~ Specify variables
fixedsamplesize = true
save = true
nsamples = 32
rng = Random.Xoshiro(42*nsamples)

if fixedsamplesize
    FILENAME = "Nfixed_" * FILENAME
    TLFILENAME = "Nfixed_" * TLFILENAME
    N = 16_000 .* ones(Int, nsamples)
else
    FILENAME = "Nvariable_" * FILENAME
    TLFILENAME = "Nvariable_" * TLFILENAME
    #~ Log-normal
    D = 10_000
    σ = 0.1
    μ = log(D) - σ^2 / 2
    Npdf = LogNormal(μ,σ)
    #~ Gamma
    # θ  = 1.0
    # Nm = 1e4 
    # α  = Nm / θ
    # Npdf = Gamma(α, θ)
    N = trunc.(Int, rand(rng, Npdf, nsamples))
end

#/ Sample from bag-of-words, both with and without replacement
dfnorep = Meris.WordSampler.samplebagofwords(N; rng=rng, replace=false)
dfrep   = Meris.WordSampler.samplebagofwords(N; rng=rng, replace=true)

#/ Sample using multinomial and multivariate hypergeometric
dfmult  = Meris.WordSampler.samplefrequencyofwords(N; rng=rng, replace=true)
dfmvhyp = Meris.WordSampler.samplefrequencyofwords(N; rng=rng, replace=false)

#/ Compute AFDs for all
znorep = Meris.AFD.compute(dfnorep)
zrep   = Meris.AFD.compute(dfrep)
zmult  = Meris.AFD.compute(dfmult)
zmvhyp = Meris.AFD.compute(dfmvhyp)
#/ Additionally compute mean and variance of the relative frequencies `x`
tlnorep = Meris.Taylor.compute(dfnorep)
tlrep   = Meris.Taylor.compute(dfrep)
tlmult  = Meris.Taylor.compute(dfmult)
tlmvhyp = Meris.Taylor.compute(dfmvhyp)

#/ Store dataframes
if save
    CSV.write(DATADIR*"noreplace_"*FILENAME, znorep)
    CSV.write(DATADIR*"replace_"*FILENAME, zrep)
    CSV.write(DATADIR*"multinomial_"*FILENAME, zmult)
    CSV.write(DATADIR*"mvhypgeom_"*FILENAME, zmvhyp)

    CSV.write(DATADIR*"noreplace_"*TLFILENAME, tlnorep)
    CSV.write(DATADIR*"replace_"*TLFILENAME, tlrep)
    CSV.write(DATADIR*"multinomial_"*TLFILENAME, tlmult)
    CSV.write(DATADIR*"mvhypgeom_"*TLFILENAME, tlmvhyp)
end
