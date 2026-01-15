#= Module to plot Heap's law for some processes and/or data =#
#/ Start module
module TaylorPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings
using FileIO, ImageTransformations

using StatsBase
using JLD2

#/ Modules
using Meris

const ICONDIR = normpath(joinpath(@__DIR__, "..", "icons/"))

#################
### FUNCTIONS ###
function plot_taylor(;
    TLDIR=Meris.DATADIR * "macro/taylor/",
    savefig=false,
    figname="tl.png"
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = 1.9 * 246
    height = width
    fig = Figure(; size=(width, height / 2), figure_padding=(2, 4, 2, 14))

    #/ Plot
    xtl = -10:0.1:10.0
    icons = []
    axes = []

    #/ RFC
    axrfc = Axis(
        fig[1, 3],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-3, 4, -3, 5)
    )
    rfc_df = JLD2.load(TLDIR * "rfc.jld2")["tldf"]
    fitrfc = fit(rfc_df)
    m, s = log10.(rfc_df[!, :omeanfrequency]), log10.(rfc_df[!, :ovarfrequency])
    m .-= mean(m)
    s .-= mean(s)
    # s .-= fitrfc.a
    rfctl = scatter!(
        axrfc, m, s,
        color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4, label=L"\text{RFCs}"
    )
    lines!(axrfc, xtl .- minimum(m) .- 0.5, 2 .* (xtl .- minimum(m)), linewidth=1, color=:black, linestyle=(:dash,:dense))
    lines!(axrfc, xtl .- minimum(m) .- 0.5, xtl, linewidth=1, color=:grey, linestyle=(:dash,:dense))
    # axislegend(axrfc, position=:rb)
    push!(icons, ICONDIR*"documents.png")
    push!(axes, (; ax=axrfc, pos=[1,3]))

    # #/ OTU
    axotu = Axis(
        fig[2, 3],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-2, 4, -4, 6)
    )
    otu_df = JLD2.load(TLDIR * "otu.jld2")["tldf"]
    fitotu = fit(otu_df)
    m, s = log10.(otu_df[!, :omeanfrequency]), log10.(otu_df[!, :ovarfrequency])
    m .-= mean(m)
    s .-= mean(s)
    # s .-= fitotu.a
    otutl = scatter!(
        axotu, m, s,
        color=:white, strokecolor=colors[2], markersize=4, strokewidth=0.4, label=L"\text{OTUs}"
    )
    lines!(axotu, xtl .- minimum(m) .- 0.5, 2 .* (xtl .- minimum(m)) , linewidth=1, color=:black, linestyle=(:dash,:dense))
    lines!(axotu, xtl .- minimum(m) .- 0.5, xtl, linewidth=1, color=:grey, linestyle=(:dash,:dense))
    # axislegend(axotu, position=:rb)
    push!(icons, ICONDIR*"bacteria.png")
    push!(axes, (; ax=axotu, pos=[2,3]))

    # #/ Lego
    axlego = Axis(
        fig[1, 4],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-1, 2, -2, 3)
    )
    lego_df = JLD2.load(TLDIR * "lego.jld2")["tldf"]
    fitlego = fit(lego_df)
    m, s = log10.(lego_df[!, :omeanfrequency]), log10.(lego_df[!, :ovarfrequency])
    m .-= mean(m)
    s .-= mean(s)
    # s .-= fitlego.a
    legotl = scatter!(
        axlego, m, s,
        color=:white, strokecolor=colors[3], markersize=4, strokewidth=0.4, label=L"\text{LEGO}"
    )
    lines!(axlego, xtl .- minimum(m) .- 0.3, 2 .* (xtl .- minimum(m)) , linewidth=1, color=:black, linestyle=(:dash,:dense))
    lines!(axlego, xtl .- minimum(m) .- 0.3, xtl, linewidth=1, color=:grey, linestyle=(:dash,:dense))
    # axislegend(axlego, position=:rb)
    push!(icons, ICONDIR*"lego.png")
    push!(axes, (; ax=axlego, pos=[1,4]))

    # #/ GTEx
    axgtex = Axis(
        fig[2, 4],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-6, 8, -10, 12)
    )
    gtex_df = JLD2.load(TLDIR * "gtex.jld2")["tldf"]
    fitgtex = fit(gtex_df)
    m, s = log10.(gtex_df[!, :omeanfrequency]), log10.(gtex_df[!, :ovarfrequency])
    m .-= mean(m)
    s .-= mean(s)
    # s .-= fitgtex.a
    gtextl = scatter!(
        axgtex, m, s,
        color=:white, strokecolor=colors[4], markersize=4, strokewidth=0.4, label=L"\text{GTEx}"
    )
    lines!(axgtex, xtl .- minimum(m) .- 1.3, 2 .* (xtl .- minimum(m)) , linewidth=1, color=:black, linestyle=(:dash,:dense))
    lines!(axgtex, xtl .- minimum(m) .- 1.3, xtl, linewidth=1, color=:grey, linestyle=(:dash,:dense))
    # axislegend(axgtex, position=:rb)
    push!(icons, ICONDIR*"gene.png")
    push!(axes, (; ax=axgtex, pos=[2,4]))

    # Plot high occupancy TL
    axtl = Axis(
        fig[1:2, 1:2],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-2, 2, -4, 4)
    )

    for (i, ax) in enumerate(axes)
        (xmin, xmax, ymin, ymax) = ax.ax.limits[]
        # lines!(ax.ax, [xmin, xmax], [ymin, ymax], linestyle=(:dash,:dense), color=:gray)
        axicon = Axis(
            fig[ax.pos...],
            width=Relative(0.22), height=Relative(0.22),
            halign=0.072, valign=0.95
        )
        icon = FileIO.load(icons[i])
        icon_small = imresize(icon, (256, 256))
        image!(axicon, rotr90(icon))
        hidedecorations!(axicon)
        hidespines!(axicon)
    end

    m, s = log10.(gtex_df[gtex_df.occupancy .> 0.99999, :][!, :omeanfrequency]), log10.(gtex_df[gtex_df.occupancy .> 0.99999, :][!, :ovarfrequency])
    m .-= mean(m)
    s .-= mean(s)
    gtextl = scatter!(
        axtl, m[1:400], s[1:400],
        color=:white, strokecolor=colors[4], markersize=4, strokewidth=0.4, label=L"\text{GTEx}"
    )

    m, s = log10.(otu_df[otu_df.occupancy .> 0.8, :][!, :omeanfrequency]), log10.(otu_df[otu_df.occupancy .> 0.8, :][!, :ovarfrequency])
    m .-= mean(m)
    s .-= mean(s)
    otutl = scatter!(
        axtl, m, s,
        color=:white, strokecolor=colors[2], markersize=4, strokewidth=0.4, label=L"\text{OTUs}"
    )
    
    m, s = log10.(rfc_df[rfc_df.occupancy .> 0.8, :][!, :omeanfrequency]), log10.(rfc_df[rfc_df.occupancy .> 0.8, :][!, :ovarfrequency])
    m .-= mean(m)
    s .-= mean(s)
    rfctl = scatter!(
        axtl, m, s,
        color=:white, strokecolor=colors[1], markersize=4, strokewidth=0.4, label=L"\text{RFCs}"
    )

    m, s = log10.(lego_df[lego_df.occupancy .> 0.1, :][!, :omeanfrequency]), log10.(lego_df[lego_df.occupancy .> 0.1, :][!, :ovarfrequency])
    m .-= mean(m)
    s .-= mean(s)
    legotl = scatter!(
        axtl, m, s,
        color=:white, strokecolor=colors[3], markersize=4, strokewidth=0.4, label=L"\text{LEGO}"
    )

    lines!(axtl, xtl, 2 .* xtl, linewidth=1, color=:black, linestyle=(:dash,:dense))
    axislegend(axtl, position=:rb, visibleframe=true)

    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

### HELPER ###
function fit(tldf)
    m = tldf[!,:meanfrequency]
    s = tldf[!,:varfrequency]
    x = log.(m)
    y = log.(s)
     #/ Fix a straight line using York's method
    #/ Calculate weights using the errors
    #~ Calculate how much of the total variation comes from presence-absence
    #  recall (σ′)² ← o⋅[σ²+μ²(1-o)], and so the ratio R = (σ′)² / (o⋅σ²) 
    o = tldf[!,:occupancy] .* (1 .- tldf[!,:occupancy])
    R = o .* tldf[!,:meanfrequency].^2 ./ tldf[!,:ovarfrequency]
    #~ filter those with ratio 0 [occupancy 0]
    sidxs = findall(x -> x > 0, R)
    x = x[sidxs]
    y = y[sidxs]
    #~ extract the errors on the mean and variance [see `taylor.jl`]
    #! note: use the δ-method to get the error on the log-transformed variables
    σx = m[sidxs] ./ m[sidxs].^2
    σy = s[sidxs] ./ s[sidxs].^2
    logcov = tldf[!,:errorcov][sidxs] ./ (m[sidxs] .* s[sidxs])
    logρ = logcov ./ sqrt.(σx .* σy)
    #~ specify the weights
    #! note: As for the line fitting only the relative weights are relevant, one could in
    #        principle scale the weights such that they are numerically more 'stable'. Yet,
    #        this may distort the error on the slope and intercept, as these are now
    #        'artificially' inflated by the weights. To bring them into a reasonable scale,
    #        we here specify the scale specifically, such that errors are reflecting the
    #        actual scatter of the means and variances and not the artificial weights.
    wx = 1.0 ./ σx
    wy = 1.0 .* sqrt.(1.0 .- R[sidxs]) ./ σy
    wscale = length(sidxs) / sum(wy)
    wx = wx .* wscale
    wy = wy .* wscale
    #~ Fit
    straightlinefit = Meris.StraightLine.weightedyorkfit(x, y, wx, wy, ρ=logρ)
    return straightlinefit
end

end # module AFDPlotter
#/ End module
