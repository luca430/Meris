#= Investigate the occupancy =#
#/ Start module
module OccuPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta
using StatsBase
using Distributions
using FHist

#/ Modules
import Meris.DATADIR as DATADIR
import Meris.StraightLine as SL

#################
### FUNCTIONS ###
function plot_occupancy(;
    DIRECTORIES = [ "rfc/", "lego/"],
    LABELS = [L"\textrm{RFC}", L"\textrm{LEGO}", L"\textrm{BCI}"],
    nbins::Int = 27,
    DIR = DATADIR * "macro/afd/",
    BASEFILENAME = "tl-stats.csv",
    savefig = false,
    figname = true
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    width = .95 * 246
    height = 3*width / 4
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"\textrm{occupancy}\;o_i", xlabelsize=12,
        ylabel=L"\textrm{log}\;R", ylabelsize=12,
        limits=(0,1,-8,1)
    )

    #/ Plot    
    for i in eachindex(DIRECTORIES)
        #/ Occupancy
        tldf = CSV.read(DIR*DIRECTORIES[i]*BASEFILENAME, DataFrame)
        #~ Rescale using the occupancy
        tldf = @chain tldf begin
            @transform(:omeanfrequency = :meanfrequency .* :occupancy)
            @transform(:ovarfrequency = :varfrequency .+ :meanfrequency.^2 .* (1 .- :occupancy))
            @transform(:ovarfrequency = :ovarfrequency .* :occupancy)
        end
        # CV = tldf[!,:varfrequency] ./ tldf[!,:meanfrequency].^2
        x = tldf[!,:occupancy]
        o = tldf[!,:occupancy] .* (1 .- tldf[!,:occupancy])
        y = log.(o .* tldf[!,:meanfrequency].^2 ./ tldf[!,:ovarfrequency])
        # x = log.(CV)
        # y = log.(tldf[!,:occupancy])
        # straightlinefit = SL.weightedyorkfit(x, y, wx, wy, ρ=tldf[!,:errorcorr])
        # @info "fit" DIRECTORIES[i] straightlinefit
        # #~ Reshuffle before plotting so it goes through (0,0)
        # l = lines!(
        #     axtl, xs, straightlinefit.b .* xs,
        #     linestyle=(:dash,:dense), linewidth=.8, color=colors[i]
        # )
        #~ Scatter w.r.t. to their weight
        # minwidth = .33
        # maxwidth = 1.0
        # strokewidth = (maxwidth - minwidth) .* wx./sum(wx) .+ minwidth
        # markersize = 4 .- (3.5 .* tldf[!,:varfrequency] ./ maximum(tldf[!,:varfrequency]))
        markersize = 2.0 .* tldf[!,:varfrequency] ./ tldf[!,:varfrequency]
        strokewidth = .6
        #~ rescale so that they have slope 1
        # yrescaled = y .- straightlinefit.a
        stl = scatter!(ax, x, y, markersize=markersize, strokewidth=strokewidth, label=LABELS[i])
        # _mean, _var = tldf[!,:meanfrequency].^2, tldf[!,:varfrequency]
        # stl = scatter!(axtl, log.(_mean), log.(_var), markersize=4, strokewidth=.7)
    end

    
    axislegend(
        ax,
        position=:rb, labelsize=9, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end


end # module OccuPlotter
#/ End module
