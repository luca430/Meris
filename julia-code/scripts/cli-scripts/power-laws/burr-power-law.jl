#/ Simple script to sample from Burr distribution and fit Pareto to tail
using FHist
using Random
using StatsBase
using JLD2

using Meris

nsamples = 10^4
rng = Random.Xoshiro(42 * nsamples)

c = 1.5
α = 1.0
λ = 2.0

logxmin = -3.
logxmax = 4.

#~ Sample and make histogram
println("Sampling and making histogram...")
# xs = Meris.Powerlaw.samplepareto(nsamples; γ=1+α*c, xmin=0.05, rng=rng)
xs = Meris.Powerlaw.sampleburr(nsamples; c=c, α=α, λ=λ, rng=rng)
zs = log10.(xs)
bins = range(logxmin, logxmax, 27)
fh = FHist.Hist1D(zs, binedges=bins) |> FHist.normalize

#~ Compute theoretical pdf
xburr = exp10.(range(logxmin,logxmax,256)) |> collect
F = Meris.ParetoLike.Burrcdf(c, α, λ)
Z = F(exp10(logxmax)) - F(exp10(logxmin))
yburr = Meris.ParetoLike.Burrpdf(xburr, c, α, λ) ./ Z

#~ Fit power law (pure Pareto) on the tail
println("Fitting pure Pareto power law to the tail...")
__zmin = median(zs)
zsorted = sort(zs)
possiblexmins = collect(range(__zmin, zsorted[end-256], 256))
paretofit = Meris.Powerlaw.fitPareto(xs, xmins=exp10.(possiblexmins))

#~ Fit generalized Pareto on the tail
println("Fitting generalized Pareto distribution to the tail...")
genparetofit = Meris.Powerlaw.fitGeneralizedPareto(xs, xmins=exp10.(possiblexmins))

#~ Compute the Hill's estimators for ξ from the data
println("Computing Hill's and log-variance estimators for ξ...")
xs_sorted = sort(xs, rev=true)
ξHill = Meris.Powerlaw.hills_estimator(xs_sorted, sorted=true)
#~ Compute the log-variance estimator for ξ from the data
ξLV = Meris.Powerlaw.log_variances(xs_sorted, sorted=true)

#~ Compute another histogram that considers only the tail
tailbins = range(log10(paretofit.xmin), logxmax, 19)
zstail = zs[zs .> log10(paretofit.xmin)]
fhtail = FHist.Hist1D(zstail, binedges=tailbins) |> FHist.normalize

#/ SAVE
println("Saving...")
OUTDIR = Meris.DATADIR * "heavy-tails/"
mkpath(OUTDIR)
OUTFILE = OUTDIR * "burr-power-law.jld2"
jldsave(
    OUTFILE;
    params=(; c=c, α=α, λ=λ),
    xs=xs,
    fh=(; x=bincenters(fh), y=fh.bincounts),
    fhtail=(; x=bincenters(fhtail), y=fhtail.bincounts),
    pdf=(; x=xburr, y=yburr),
    paretofit=paretofit,
    genparetofit=genparetofit,
    ξHill = ξHill,
    ξLV = ξLV
)
