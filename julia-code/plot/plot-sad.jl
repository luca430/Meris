#= Module to plot SADs/CADs for some system parameters =#
#/ Start module
module SADPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings
# using LsqFit
using JLD2
using StatsBase

#/ Local packages, modules, and directories
import Meris.DATADIR as DATADIR

#################
### FUNCTIONS ###
function plot_synthetic_sad(;
    γv = [0.5, 1.5, 2.5],
    DIR=DATADIR*"sad/synthetic/",
    savefig=false,
    figname=nothing
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(width,height), figure_padding=(2,8,2,8))
    ax = Axis(
        fig[1,1],
        xlabel=L"\textrm{aggregated frequencies}\;\log_{10}\tilde{\nu}", xlabelsize=11,
        ylabel=L"\textrm{density}\;\log_{10}\,p(\tilde{\nu})", ylabelsize=11,
        limits=(-7,-2,-6,1)
    )

    for i in eachindex(γv)
        filename = "synthetic-sad-gamma$(γv[i]).jld2"
        result = JLD2.load(DIR*filename)
        
        scatter!(
            ax, result["n"], log10.(result["p"]), markersize=4, color=:white,
            strokecolor=(colors[i],0.6), strokewidth=.42
        )

        lines!(
            ax, log10.(result["paretox"]), log10.(result["paretoy"]), linestyle=:dash,
            linewidth=.8, color=colors[i], label=L"\gamma=%$(γv[i])"
        )
    end

    #/ Legend
    axislegend(
        ax,
        position=:rt, labelsize=9, patchsize=(8,20),
        margin=(2,8,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module SADPlotter

#/ End module
