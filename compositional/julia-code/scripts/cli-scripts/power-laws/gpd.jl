#/ Simple script to sample from generalized Pareto distribution and fit generalized Pareto to tail
using Distributions
using FHist
using Random
using StatsBase
using JLD2

using Meris

nsamples = 10^4
rng = Random.Xoshiro(42 * nsamples)

μ = 1.0
σ = 1.5
ξ = 1.667

logxmin = 0.
logxmax = 5.

#~ Sample and make histogram
println("Sampling and making histogram...")
# xs = Meris.Powerlaw.samplepareto(nsamples; γ=1+α*c, xmin=0.05, rng=rng)
gpd = Distributions.GeneralizedPareto(μ, σ, ξ)
xs = rand(rng, gpd, nsamples)
zs = log10.(xs)
bins = range(logxmin, logxmax, 27)
fh = FHist.Hist1D(zs, binedges=bins) |> FHist.normalize

#~ Compute theoretical pdf
xgpd = exp10.(range(logxmin,logxmax,256)) |> collect
F = Meris.ParetoLike.generalizedParetocdf(σ, ξ; xmin=μ)
Z = F(exp10(logxmax)) - F(exp10(logxmin))
ygpd = Meris.ParetoLike.generalizedParetopdf(xgpd, σ, ξ; xmin=μ) ./ Z

#~ Fit power law on the tail
println("Fitting pure Pareto power law to the tail...")
__zmin = median(zs)
zsorted = sort(zs)
possiblexmins = collect(range(__zmin, zsorted[end-256], 256))
paretofit = Meris.Powerlaw.fitPareto(xs, xmins=exp10.(possiblexmins))

println("Fitting generalized Pareto distribution to the tail...")
genparetofit = Meris.Powerlaw.fitGeneralizedPareto(xs, xmins=exp10.(possiblexmins))
eygpd = Meris.ParetoLike.generalizedParetopdf(
    xgpd, genparetofit.σ, genparetofit.ξ; xmin=genparetofit.xmin
)

#~ Compute another histogram that considers only the tail
tailbins = range(log10(paretofit.xmin), logxmax, 19)
zstail = zs[zs .> log10(paretofit.xmin)]
fhtail = FHist.Hist1D(zstail, binedges=tailbins) |> FHist.normalize

#/ SAVE
println("Saving...")
OUTDIR = Meris.DATADIR * "heavy-tails/"
mkpath(OUTDIR)
OUTFILE = OUTDIR * "gpd-power-law.jld2"
jldsave(
    OUTFILE;
    params=(; μ=μ, σ=σ, ξ=ξ),
    xs=xs,
    fh=(; x=bincenters(fh), y=fh.bincounts),
    fhtail=(; x=bincenters(fhtail), y=fhtail.bincounts),
    pdf=(; x=xgpd, y=ygpd),
    epdf=(; x=xgpd, y=eygpd),
    paretofit=paretofit,
    generalizedparetofit=genparetofit
)
