#= Module to plot (distributions of) document sizes =#
#/ Start module
module SizePlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings
using CSV, DataFrames
using StatsBase
using Distributions
using FHist

#/ Modules
import Meris.DATADIR as DATADIR     #~ Data with (parsed/analyzed) results, for plotting
import Meris.ARXIVDIR as RARXIVDIR  #~ Directory with raw arXiv data

#################
### FUNCTIONS ###
function plot_documentsize(;
    LEGODIR  = DATADIR * "documentsize/lego/",
    ARXIVDIR = DATADIR * "documentsize/arxiv/",
    ARXIVCATEGORIES = readlines(RARXIVDIR*"categories.txt"),
    nbins=17
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
        xlabel=L"\textrm{document size}\;\log\,N", xlabelsize=11,
        ylabel=L"\textrm{pdf}\;P(N)", ylabelsize=11,
        yscale=log10,
        limits=(-4,4,1e-4,1)
    )

    #/ Lego sets
    legofilename = "lego-documentsize.csv"
    legodf = CSV.read(LEGODIR*legofilename, DataFrame)
    z = log.(legodf.documentsize)
    μ = mean(log.(z))
    σ = std(log.(z))
    z = (log.(z) .- μ) ./ σ
    # gamma = Distributions.Gamma(4.0, 1.0)
    # z = log.(rand(gamma, 4*256))
    zmin, zmax = extrema(z)
    binedges = range(zmin, zmax, nbins)
    fh = FHist.Hist1D(z; counttype=Int, binedges=binedges, overflow=true) |> normalize
    # return fh
    s = scatter!(ax, bincenters(fh), fh.bincounts, label=L"\textrm{LEGO}")
    #~ naively fit a Gamma distribution
    gamma = Distributions.fit_mle(Gamma, exp.(z))
    xgam = exp.(zmin:0.01:zmax)
    ygam = xgam .* Distributions.pdf.(gamma, xgam)
    lines!(ax, log.(xgam), ygam, linewidth=1.)

    


    #/ arXiv papers
    # N = Int[]
    # for CATEGORY in ARXIVCATEGORIES
    #     arxivfilename = "arxiv-$(CATEGORY)-documentsize.csv"
    #     df = CSV.read(ARXIVDIR*arxivfilename, DataFrame)
    #     df = filter(:totalcounts => n -> n .> 1_000, df)
    #     # push!(N, df.totalcounts...)
    #     z = log.(df.totalcounts)
    #     # μ = mean(log.(N))
    #     # σ = std(log.(N))
    #     # z = (log.(N) .- μ) ./ σ
    #     zmin, zmax = extrema(z)
    #     binedges = range(zmin, zmax, nbins)
    #     fh = FHist.Hist1D(z; counttype=Int, binedges=binedges, overflow=true) |> normalize
      
    #     gamma = Distributions.fit_mle(Gamma, exp.(z))
    #     xgam = exp.(zmin:0.01:zmax)
    #     ygam = xgam .* Distributions.pdf.(gamma, xgam)
    #     lines!(ax, log.(xgam), log.(ygam), linewidth=1.)
    #     #~ Scatter
    #     categorylabel = replace(CATEGORY, "-" => "‐")
    #     s = scatter!(ax, bincenters(fh), log.(fh.bincounts), label=L"\textrm{%$(categorylabel)}")
    # end
    #/ Add legend    
    axislegend(
        ax,
        position=:lt, labelsize=9, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    
	  return fig
end

end # module SizePlotter
#/ End module
