#= Simple module to test double Pareto fitting procedures on Wikitext =#

#/ Start module
module Wiki

using StatsBase
using FHist

using Meris

#################
### FUNCTIONS ###
function wikifit(; nbins=31)
    #~ Load data & compute (relative) frequencies
    words = Meris.WikitextSampler.load_bagofwords()
    counts = values(countmap(words))
    x = counts ./ sum(counts)
    
    #~ Compute the histogram
    logx = log10.(x)
    bmin, bmax = extrema(logx)
    binedges = range(bmin, bmax, nbins)
    fh = FHist.Hist1D(logx; binedges=binedges, overflow=true) |> normalize
    #~ Specify all ε's (xmin's) to try
    xmins = exp10.(range(bmin, bmax, 64))
    
    #~ Fit some pure Pareto on the right-tail
    pfit = Meris.Powerlaw.fitPareto(x; xmins=xmins)
    Pareto = pfit.Pareto
    xpdf = exp10.(range(log10(Pareto.ε), exp10(bmax), 256))
    ypdf = (10 .^ log10.(xpdf)) .* Meris.ParetoDistribution.pdf.(Pareto, xpdf)
    ztail = sum(x .> Pareto.ε) / length(x)
    ypdf = ypdf .* (log(10) * ztail)
    pareto = (; xpdf=xpdf, ypdf=ypdf)

    #~ Using this fit, fit a bounded Pareto on the left-part
    #  note: in the Wikitext example it may not make the most sense
    bx = x[x .< Pareto.ε]
    bxmins = exp10.(range(bmin, log10(Pareto.ε), 64))
    bpfit = Meris.Powerlaw.fitBoundedPareto(bx; xmins=bxmins, εmax=Pareto.ε)
    bPareto = bpfit.BoundedPareto
    xpdf = exp10.(range(log10(bPareto.ε), exp10(bPareto.εmax), 256))
    ypdf = (10 .^ log10.(xpdf)) .* Meris.ParetoDistribution.pdf.(bPareto, xpdf)
    ztail = sum(bx .> bPareto.ε) / length(bx)
    ypdf = ypdf .* (log(10) * ztail)
    bpareto = (; xpdf=xpdf, ypdf=ypdf)    
    
    #~ Fit some Pareto on the right-tail  
    gpfit = Meris.Powerlaw.fitGeneralizedPareto(x; xmins=xmins)
    gPareto = gpfit.GeneralizedPareto
    #~ Compute proper transformation of p(x) from the fit
    xpdf = exp10.(range(log10(gPareto.ε), exp10(bmax), 256))
    ypdf = (10 .^ log10.(xpdf)) .* Meris.ParetoDistribution.pdf.(gPareto, xpdf)
    ztail = sum(x .> gPareto.ε) / length(x)
    ypdf = ypdf .* (log(10) * ztail)
    gpareto = (; xpdf=xpdf, ypdf=ypdf)
    
    return (; fh=fh, pareto=pareto, bpareto=bpareto, gpareto=gpareto)
end


end # module Wiki
#/ End module
