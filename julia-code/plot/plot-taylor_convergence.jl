#= Module to plot Taylor's law exponent as a function of occupancy and sample size =#
#/ Start module
module TaylorConvPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta, JLD2
using StatsBase
using Distributions
using FHist

#/ Modules
import Meris

#################
### FUNCTIONS ###
function plot(;
        datadir = Meris.DATADIR * "macro/",
        savefig = false,
        figname = true
    )

    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(2*width,height), figure_padding=(2,4,2,14))
    
    ax1 = Axis(
        fig[1,1],
        xlabel=L"\text{occupancy } o",
        ylabel=L"b"
    )
    
    ax2 = Axis(
        fig[1,2],
        xlabel=L"\text{sample size } N",
        ylabel=L"b",
        xscale=log10
    )

    @load datadir * "b_o_dict_genetic.jld2" d
    for (k, v) in d
        scatter!(ax1, v[1], v[2], markersize=5, label=L"\text{%$k}")
        lines!(ax1, v[1], v[2], linewidth=0.8)
    end

    @load datadir * "b_N_dict_genetic.jld2" d
    for (k, v) in d
        scatter!(ax2, v[1], v[2], markersize=5)
        lines!(ax2, v[1], v[2], linewidth=0.8)
    end

    leg = Legend(fig[1, 3], ax1; tellheight=false)

    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))

    return fig
end


########################
### HELPER FUNCTIONS ###


##############################

end # module TaylorConvPlotter
#/ End module
