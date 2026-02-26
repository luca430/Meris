#= Simple script to compute the survival function for a bag-of-legos =#
using JLD2
using DataFrames, DataFramesMeta
using FHist

using Meris
DATADIR = Meris.DATADIR * "survival/lego/"
mkpath(DATADIR)

#~ Specify variables
save = true
minquantity = 100
mindistinctpieces = 30
nbins = 31

#~ Get the bag-of-legos and compute relative frequencies
legodf = Meris.LegoSampler.filterlegos(;
    minquantity=minquantity, mindistinctpieces=mindistinctpieces
)
@transform!(legodf, :frequency = :counts ./ :nreads)
bagoflegos = sort(@select(legodf, :component_id, :frequency), :frequency)

#/ Compute histogram of log-frequencies
println("Computing histogram...")
logfreq = log.(bagoflegos[!,:frequency])
fmin, fmax = extrema(logfreq)
edges = collect(range(fmin, fmax, nbins))
fh = FHist.Hist1D(logfreq; binedges=edges, overflow=true) |> normalize

#/ Compute survival function S(t) = Pr[x > t]
println("Computing survival function...")
n = nrow(bagoflegos)
unqfreqs = unique(bagoflegos[!,:frequency])
#~ For each unique frequency ν, get the last index k-1 for which x < ν.
#  Then the fraction of elements larger than ν is exactly (n-k)/n
k = searchsortedlast.(Ref(bagoflegos[!,:frequency]), unqfreqs) .+ 1
S = (n .- k .+ 1) ./ n

#/ Do quick-n-dirty hypothesis check to check for Pareto
println("Fitting Pareto and checking untruncated Pareto hypothesis...")
X = unqfreqs[end]
possiblexmins = exp10.(range(fmin, fmax-1, 256))
paretofit = Meris.Powerlaw.fitPareto(bagoflegos[!,:frequency], xmins=possiblexmins)
C = paretofit.xmin^paretofit.γ
qs = [0.01,0.05,0.1,0.2]
ps = similar(qs, Float64)
ispareto = zeros(Bool, length(qs))
for i in eachindex(qs)
    ispareto[i] = X > (n * C / -log(qs[i]))^(1 / paretofit.γ)
    ps[i] = exp(-n*C*X^(-paretofit.γ))
end

#/ Save
if save
    filename = "lego-survival.jld2"
    jldsave(DATADIR*filename;
        fh=fh,
        isPareto=(; qs=qs, ps=ps, ispareto=ispareto),
        survivalfunction=(; t=unqfreqs, S=S)
    )
end
