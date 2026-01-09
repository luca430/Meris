#= Module to plot Zipf's law for some processes and/or data =#
#/ Start module
module ZipfPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings
using LsqFit
using JLD2
using StatsBase

#/ Modules
import Meris.DATADIR as DATADIR     #~ Data with (parsed/analyzed) results, for plotting
import Meris.ARXIVDIR as RARXIVDIR  #~ Directory with raw arXiv data

#################
### FUNCTIONS ###
function plot_zipf(;
    LEGODIR = DATADIR * "zipf/lego/",
    TREEDIR = DATADIR * "zipf/bci.tree/",
    RFCDIR  = DATADIR * "zipf/rfc/",
    OTUDIR  = DATADIR * "zipf/otu/",
    savefig = false,
    figname = true
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = .95 * 246 * 1.42
    height = 3*width / 4.67
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"\textrm{rank}\;\log_{10}r", xlabelsize=11,
        ylabel=L"\textrm{frequency}\;\log_{10} \nu", ylabelsize=11,
        # limits=(1,6,0,5)
    )
    #~ Define some functions for fitting
    logfun(x,p) = @. p[1] + p[2]*log(x)
    powfun(x,p) = @. p[1] + p[2]*x^p[3]

    #/ Load and plot LEGO
    legofilename = "lego-zipf.jld2"
    legodf = JLD2.load(LEGODIR*legofilename)
    legorank, legofreq = legodf["rank"], legodf["frequency"]
    scatter!(
        ax, log10.(legorank), log10.(legofreq),
        marker=:rect, markersize=4, color=:white, strokecolor=colors[1], strokewidth=.42,
        label=L"\textrm{LEGO}"
    )
    #/ Load and plot BCI trees
    bcitreefilename = "bcitree-zipf.jld2"
    treedf = JLD2.load(TREEDIR*bcitreefilename)
    treerank, treefreq = treedf["rank"], treedf["frequency"]
    scatter!(
        ax, log10.(treerank), log10.(treefreq),
        marker=:xcross, markersize=4, color=:white, strokecolor=colors[2], strokewidth=.42,
        label=L"\textrm{bci.tree}"
    )
    #/ Load vocabulary OTUs
    environments = ["gut1", "gut2", "seawater"]
    for (i, env) in enumerate(environments)
        otudf = JLD2.load(OTUDIR*"otu-$(env)-zipf.jld2")
        oturank, otufreq = otudf["rank"], otudf["frequency"]
        scatter!(
            ax, log10.(oturank), log10.(otufreq),
            marker=:circle, markersize=4, color=:white, strokecolor=colors[2+i], strokewidth=.42,
            label=L"\textrm{%$(env)}"
        )
    end
    
    #/ Load and plot RFC
    rfcfilename = "rfc-zipf.jld2"
    rfcdf = JLD2.load(RFCDIR*rfcfilename)
    rfcrank, rfcfreq = rfcdf["rank"], rfcdf["frequency"]
    scatter!(
        ax, log10.(rfcrank), log10.(rfcfreq),
        marker=:utriangle, markersize=4, color=:white,
        strokecolor=colors[3+length(environments)], strokewidth=.42,
        label=L"\textrm{RFC}"
    )
    
    #/ Add legend    
    axislegend(
        ax,
        position=:lb, labelsize=9, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module HeapPlotter
#/ End module
