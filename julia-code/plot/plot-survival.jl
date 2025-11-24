#= Survival plots =#
#/ Start module
module SurvivalPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings
using StatsBase
using JLD2

#/ Modules
using Meris

#################
### FUNCTIONS ###
function plot(;
    DATADIR = Meris.DATADIR * "survival/lego/",
    FILENAME = "lego-survival.jld2",
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)
    #/ Make figure
    width = .9 * 246
    height = 3*width / 4.37
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"\log_{10}(t)", xlabelsize=11,
        ylabel=L"\textrm{Pr}[x > t]", ylabelsize=11,
        limits=(-4,0,-6,1)
    )

    #/ Load data
    jldata = JLD2.load(DATADIR*FILENAME)
    Surv = jldata["survivalfunction"]
    isPareto = jldata["isPareto"]
    print(isPareto)

    lines!(ax, log10.(Surv.t), log10.(Surv.S), linewidth=1.)

    return fig
end
end # module SurvivalPlotter
#/ End module
