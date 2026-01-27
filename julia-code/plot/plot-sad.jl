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

function plot_synthetic_sadN(;
    log10N = [7],
    ids = [1, 2],
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
        xlabel=L"\textrm{rel. abundance}\;\log_{10}\nu", xlabelsize=11,
        ylabel=L"\textrm{density}\;\log_{10}\,p(\nu)", ylabelsize=11,
        limits=(-10,0,-4,1)
    )
    

    for i in eachindex(log10N), k in ids
        filename = "synthetic-sad-log10N$(log10N[i])-$(k).jld2"
        result = JLD2.load(DIR*filename)

        (i == 1) && (vlines!(ax, [-6-log10(1.5)], color=:gray, linestyle=:dot, linewidth=.5))
        
        #~ Pareto
        if i == length(log10N) & k == 1
            lines!(
                ax, log10.(result["paretox"]), log10.(result["paretoy"]), linestyle=:dash,
                linewidth=.8, color=:black, 
            )
        end
        
        scatter!(
            ax, result["n"], log10.(result["p"]), markersize=4, color=:white,
            strokecolor=colors[i], strokewidth=.42, label=L"N=10^{%$(log10N[i])}"
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

function plot_synthetic_sad_illustration(;
    log10N = 7,
    ids = [1],
    γ = 0.5,
    DIR=DATADIR*"sad/synthetic/",
    savefig=false,
    figname=nothing
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = [:black, :grey]
    fcolors = [(:grey, 0.4), (:grey, 0.6)]

    width = .45 * 246
    height = width
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"\textrm{frequency}\;\log_{10}\nu", xlabelsize=10,
        ylabel=L"\textrm{density}\;\log_{10}\,p(\nu)", ylabelsize=10,
        limits=(-8,0,-4,0),
        xminorticks=IntervalsBetween(4), yminorticks=IntervalsBetween(4)
    )
    

    for k in ids
        filename = "synthetic-sad-log10N$(log10N)-$(k).jld2"
        result = JLD2.load(DIR*filename)        
        #~ Pareto
        if k == 1
            lines!(
                ax, log10.(result["paretox"]), log10.(result["paretoy"]),
                linestyle=(:dash,:dense), linewidth=1., color=:black, 
            )
        end
        #~ Scatter
        scatter!(
            ax, result["n"], log10.(result["p"]), markersize=3, color=fcolors[k],
            strokecolor=(:black,0.8), strokewidth=.42, label=L"N_%$(k)"
        )
    end

    #/ Clarifying labels
    text!(
        -4.5,-1.8, rotation=-_get_angle(ax,γ), text=L"\propto \nu^{-\gamma}",
        align=(:left,:top), fontsize=12
    )

    #/ Legend
    # axislegend(
    #     ax,
    #     position=:rt, labelsize=10, patchsize=(8,20),
    #     margin=(2,8,0,0), patchlabelgap=2, padding=(0,0,0,0)
    # )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

########################
### HELPER FUNCTIONS ###
function _get_angle(ax, γ)
    (xmin, xmax, ymin, ymax) = ax.limits[]
    sx = 1 / (xmax - xmin)
    sy = γ / (ymax - ymin)
    angle = atan(sy, sx)
    return angle
end

end # module SADPlotter

#/ End module
