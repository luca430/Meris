#= Module to plot Heap's law for some processes and/or data =#
#/ Start module
module AFDPlotter

#/ Packages
using CairoMakie
using MakiePublication
using CSV, DataFrames
using FHist
using LaTeXStrings
using SpecialFunctions, Optim
# using LsqFit
using JLD2
using StatsBase

#/ Modules
import Meris.DATADIR as DATADIR

#################
### FUNCTIONS ###
function plot_brownafd(;
    methods = ["noreplace", "replace", "multinomial", "mvhypgeom"],
    nbins::Int = 19,
    DIR = DATADIR * "macro/afd/",
    FILENAME = "Nfixed_zvalues.csv",
    TLFILENAME = "Nfixed_stats.csv",
    savefig = false,
    figname = true
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"z", xlabelsize=11,
        ylabel=L"p(z)", ylabelsize=11,
        # yscale=log10,
        limits = (-5.,3.,0,0.6)
    )

    #/ Plot
    β = zeros(length(methods))
    for i in eachindex(methods)
        zdf = CSV.read(DIR*methods[i]*"_"*FILENAME, DataFrame)
        z = zdf[!,:z]
        zmin, zmax = extrema(zdf[!,:z])
        binedges = range(zmin, zmax, nbins)
        fh = FHist.Hist1D(z; counttype=Int, binedges=binedges, overflow=true) |> normalize
        #~ Compute β using Taylor's law
        tldf = CSV.read(DIR*methods[i]*"_"*TLFILENAME, DataFrame)
        c = mean(tldf[!,:meanfrequency] .^2 ./ tldf[!,:varfrequency])
        β[i] = 2*c/(c-1)
        #~ Try to plot log-Lomax
        #~ Find the `b` param of the log-lomax by fitting [@TODO Use Taylor's law instead]
        # neglogll(p) = -sum(loglomax(z, p))
        # res = Optim.optimize(neglogll, 1e-6, 1e4)
        # btemp = Optim.minimizer(res)
        # b += btemp
        
        #~ Scatter
        s = scatter!(ax, bincenters(fh), fh.bincounts, label=L"\textrm{%$(methods[i])}")
    end

    βplot = mean(β)
    zplot = -5:0.01:3
    pplot = loglomax.(zplot, βplot)
    lines!(ax, zplot, exp.(pplot), linewidth=.8, color=:black, label=L"\textrm{log-Lomax}")

    
    axislegend(
        ax,
        position=:lt, labelsize=10, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    
    return fig
end

########################
### HELPER FUNCTIONS ###

#/ A ZOO OF DISTRIBUTIONS
function loggamma(z, α)
    return α*sqrt(trigamma(α)) .* z .+ α*digamma(α) .- exp.(z .* sqrt(trigamma(α)) .+ digamma(α)) .+ 0.5*log(trigamma(α)) .- loggamma(α)
end

"Logarithmic Lomax distribution"
function loglomax(z, b)
    s = sqrt(trigamma(1) + trigamma(b))
    m = digamma(1) - digamma(b)
    return log(s * b) .+ z .* s .+ m .- (b + 1) .* log.(1 .+ exp.(z .* s .+ m))
end

function lrln(z, σ)
    return -z .^ 2 ./ 2 .- log(sqrt(σ^2 * 2 * π))
end
##############################

end # module AFDPlotter
#/ End module
