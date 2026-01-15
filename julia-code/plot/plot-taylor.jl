#= Module to plot Heap's law for some processes and/or data =#
#/ Start module
module TLPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using StatsBase
using JLD2

#/ Modules
import Meris.DATADIR as DATADIR
import Meris.StraightLine as SL

#################
### FUNCTIONS ###
function plot_taylor(;
    TLDIR=DATADIR * "taylor/",
    rescale=false,
    savefig=false,
    figname=true
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = 1.9 * 246
    height = width
    fig = Figure(; size=(width, height / 2), figure_padding=(2, 4, 2, 14))

    #/ Plot
    # xtl = -20:0.1:0.0
    # ytl = copy(xtl).*2
    # lines!(ax, xtl, ytl, linewidth=1, color=:black, linestyle=(:solid,:dense))

    #/ Gutemberg
    axtl = Axis(fig[1:2, 1:2])
    axrfc = Axis(
        fig[1, 3],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-6, -1, -10, -2)
    )
    gutemberg_df = JLD2.load(TLDIR * "")
    m, s = gutemberg_df["omean"], gutemberg_df["ovar"]
    gutemberg_df = scatter!(
        axrfc, log10.(m), log10.(s),
        color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4
    )

    #/ OTU
    axotu = Axis(
        fig[2, 3],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-6, -1, -10, -2)
    )
    otu_df = JLD2.load(TLDIR * "")
    m, s = otu_df["omean"], otu_df["ovar"]
    otutl = scatter!(
        axotu, log10.(m), log10.(s),
        color=:white, strokecolor=colors[2], markersize=4, strokewidth=0.4
    )

    #/ Lego
    axlego = Axis(
        fig[1, 4],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-6, -1, -10, -2)
    )
    lego_df = JLD2.load(TLDIR * "")
    m, s = lego_df["omean"], lego_df["ovar"]
    legotl = scatter!(
        axlego, log10.(m), log10.(s),
        color=:white, strokecolor=colors[3], markersize=4, strokewidth=0.4
    )

    #/ GTEx
    axgtex = Axis(
        fig[2, 4],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-6, -1, -10, -2)
    )
    gtex_df = JLD2.load(TLDIR * "")
    m, s = gtex_df["omean"], gtex_df["ovar"]
    gtextl = scatter!(
        axgtex, log10.(m), log10.(s),
        color=:white, strokecolor=colors[4], markersize=4, strokewidth=0.4
    )

    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module AFDPlotter
#/ End module
