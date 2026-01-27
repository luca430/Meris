#= Module to plot Heap's law for some processes and/or data =#
#/ Start module
module TaylorPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using FileIO, ImageTransformations
using StatsBase, Random
using JLD2

#/ Modules
import Meris.DATADIR as DATADIR
import Meris.FIGUREDIR as FIGDIR
import Meris.StraightLine as SL

#################
### FUNCTIONS ###
function plot_taylor(;
    TLDIR = DATADIR * "taylor/",
    ICONDIR = FIGDIR * "icons/",
    rescale = false,
    savefig = false,
    figname = true
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = 1.9 * 246
    height = width
    fig = Figure(; size=(width,height/2), figure_padding=(2,4,2,14))

    icons = []
    axes = []
    
    #/ [...]
    axtl = Axis(fig[1:2,1:2])
    
    #/ Plot Taylor's law for RFC
    axrfc = Axis(
        fig[1,3],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-8,-1,-14,0)
    )
    rfcdf = JLD2.load(TLDIR*"rfc/rfc-taylor.jld2")
    m, s = rfcdf["omean"], rfcdf["ovar"]
    rfctl = scatter!(
        axrfc, log10.(m), log10.(s),
        color=:white, strokecolor=colors[1], markersize=4, strokewidth=.4
    )
    # push!(labels, L"\textrm{RFC}")
    push!(icons, ICONDIR*"documents.png")
    push!(axes, (; ax=axrfc, pos=[1,3]))
    
    #/ OTU
    axotu = Axis(
        fig[1,4],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-5,-1,-10,-2)
    )
    otudf = JLD2.load(TLDIR*"otu/otu-gut1-taylor.jld2")
    m, s = otudf["omean"], otudf["ovar"]
    otutl = scatter!(
        axotu, log10.(m), log10.(s),
        color=:white, strokecolor=colors[2], markersize=4, strokewidth=.4
    )
    # push!(labels, L"\textrm{OTU}")
    push!(icons, ICONDIR*"bacteria.png")
    push!(axes, (; ax=axotu, pos=[1,4]))
    # return axotu

    #/ Lego
    axlego = Axis(
        fig[2,3],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-8,-2,-12,0)
    )
    legodf = JLD2.load(TLDIR*"lego/lego-taylor.jld2")
    m, s = legodf["omean"], legodf["ovar"]
    legotl = scatter!(
        axlego, log10.(m), log10.(s),
        color=:white, strokecolor=colors[3], markersize=4, strokewidth=.4
    )
    # push!(labels, L"\textrm{Lego}")
    push!(icons, ICONDIR*"lego.png")
    push!(axes, (; ax=axlego, pos=[2,3]))

    #/ Genes
    axgtex = Axis(
        fig[2,4],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=11,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=11,
        limits=(-10,0,-20,0)
    )
    gtexdf = JLD2.load(TLDIR*"gtex/gtex-taylor.jld2")
    m, s = gtexdf["omean"], gtexdf["ovar"]
    gtextl = scatter!(
        axgtex, log10.(m), log10.(s),
        color=:white, strokecolor=colors[4], markersize=4, strokewidth=.4
    )
    push!(icons, ICONDIR*"gene.png")
    push!(axes, (; ax=axgtex, pos=[2,4]))

    #/ Plot lines and icons
    for (i, ax) in enumerate(axes)
        (xmin, xmax, ymin, ymax) = ax.ax.limits[]
        lines!(ax.ax, [xmin, xmax], [ymin, ymax], linestyle=(:dash,:dense), color=:gray)
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
    
    return fig
    
    
    # for i in eachindex(DIRECTORIES)
    #     #/ Taylor's law
    #     tldf = CSV.read(DIR*DIRECTORIES[i]*BASETLFILENAME, DataFrame)
    #     m = tldf[!,:meanfrequency]
    #     s = tldf[!,:varfrequency]
    #     x = log.(m)
    #     y = log.(s)

    #     #~ Compute rescaled moments using the occupancy        
    #     tldf = @chain tldf begin 
    #         @transform(:omeanfrequency = :meanfrequency .* :occupancy)
    #         @transform(:ovarfrequency = :varfrequency .+ (1 .- :occupancy) .* :meanfrequency.^2)
    #         @transform(:ovarfrequency = :ovarfrequency .* :occupancy)
    #     end
    #     #~ Use them if `rescale=true`
    #     if rescale
    #         x = log.(tldf[!,:omeanfrequency])
    #         y = log.(tldf[!,:ovarfrequency])
    #     end
        
    #     #/ Fix a straight line using York's method
    #     #/ Calculate weights using the errors
    #     #~ Calculate how much of the total variation comes from presence-absence
    #     #  recall (σ′)² ← o⋅[σ²+μ²(1-o)], and so the ratio R = (σ′)² / (o⋅σ²) 
    #     o = tldf[!,:occupancy] .* (1 .- tldf[!,:occupancy])
    #     R = o .* tldf[!,:meanfrequency].^2 ./ tldf[!,:ovarfrequency]
    #     #~ filter those with ratio 0 [occupancy 0]
    #     sidxs = findall(x -> x > 0, R)
    #     x = x[sidxs]
    #     y = y[sidxs]
    #     #~ extract the errors on the mean and variance [see `taylor.jl`]
    #     #! note: use the δ-method to get the error on the log-transformed variables
    #     σx = m[sidxs] ./ m[sidxs].^2
    #     σy = s[sidxs] ./ s[sidxs].^2
    #     logcov = tldf[!,:errorcov][sidxs] ./ (m[sidxs] .* s[sidxs])
    #     logρ = logcov ./ sqrt.(σx .* σy)
    #     #~ specify the weights
    #     #! note: As for the line fitting only the relative weights are relevant, one could in
    #     #        principle scale the weights such that they are numerically more 'stable'. Yet,
    #     #        this may distort the error on the slope and intercept, as these are now
    #     #        'artificially' inflated by the weights. To bring them into a reasonable scale,
    #     #        we here specify the scale specifically, such that errors are reflecting the
    #     #        actual scatter of the means and variances and not the artificial weights.
    #     wx = 1.0 ./ σx
    #     wy = 1.0 .* sqrt.(1.0 .- R[sidxs]) ./ σy
    #     wscale = length(sidxs) / sum(wy)
    #     wx = wx .* wscale
    #     wy = wy .* wscale
    #     #~ Fit
    #     straightlinefit = SL.weightedyorkfit(x, y, wx, wy, ρ=logρ)
    #     @info "fit" DIRECTORIES[i] straightlinefit
    #     #~ Reshuffle before plotting so it goes through (0,0)        
    #     xs = -20:1.0:0.0
    #     bs = round(straightlinefit.b, sigdigits=3)
    #     σs = round(straightlinefit.σb, sigdigits=2)
    #     blabel = L"b = %$(bs)\,(%$(σs))"
    #     l = lines!(
    #         axtl, xs, straightlinefit.b .* xs, label=blabel,
    #         linestyle=(:dash,:dense), linewidth=.8, color=colors[i]
    #     )
    #     #~ Scatter w.r.t. to their weight
    #     minwidth = .33
    #     maxwidth = 1.0
    #     strokewidth = (maxwidth - minwidth) .* wx./sum(wx) .+ minwidth
    #     #~ shift so that they go through the origin
    #     yshifted = y .- straightlinefit.a
    #     stl = scatter!(axtl, x, yshifted, markersize=4, strokewidth=strokewidth)
    #     # _mean, _var = tldf[!,:meanfrequency].^2, tldf[!,:varfrequency]
    #     # stl = scatter!(axtl, log.(_mean), log.(_var), markersize=4, strokewidth=.7)
    # end

    
    # axislegend(
    #     ax,
    #     position=:lt, labelsize=9, patchsize=(8,20),
    #     margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    # )
    # axislegend(
    #     axtl,
    #     position=:lt, labelsize=9, patchsize=(8,20),
    #     margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    # )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

function plot_synthetictaylor(;
    TLDIR = DATADIR * "taylor/synthetic/",
    filename = "synthetic-taylor.jld2",
    savefig = false,
    figname = nothing,
    n = 42,
    rng = Random.Xoshiro(42)
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)
    colors = MakiePublication.COLORS[begin]

    width = .45 * 246
    height = width
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    
    #/ Plot Taylor's law for synthetic data
    ax = Axis(
        fig[1,1], aspect=1,
        xlabel=L"\textrm{mean}\;\log_{10}\,\mu", xlabelsize=11,
        ylabel=L"\textrm{variance}\;\log_{10}\,\sigma^2", ylabelsize=11,
        xminorticks=IntervalsBetween(4),
        yminorticks=IntervalsBetween(4),
        limits=(-4,4,-4,6)
    )

    (xmin, xmax, ymin, ymax) = ax.limits[]
    xband = [0, 4]
    band!(
        ax, xband, ymin .* ones(length(xband)), ymax .* ones(length(xband)),
        color=(:gray, 0.2)
    )

    #~ Load data
    db = JLD2.load(TLDIR*filename)
    logm = log10.(filter(x->x>0, db["mean"]))
    logs = log10.(filter(x->x>0, db["var"]))
    #~ Take a [small] subsample
    nsamples = length(logm)
    _order = sortperm(logm)
    logm = logm[_order]
    logs = logs[_order]
    mmin, mmax = extrema(logm)
    idxs = Array{Int}(undef, n)
    for (i, m) in enumerate(range(mmin, mmax, n))
        (m == mmin) && (idxs[begin] = 1; continue)
        (m == mmax) && (idxs[end] = nsamples; break)
        idxs[i] = idxs[i-1] + findfirst(x -> x > m, logm[idxs[i-1]:end])
    end
    mplot = logm[idxs]
    splot = logs[idxs]
    
    #/ Scatter synthetic data
    scatter!(
        ax, mplot, splot,
        color=(colors[1],0.7), strokecolor=:black, markersize=4, strokewidth=.3
    )
    #/ Lines    
    #~ Determine and plot line with b=1
    a1 = minimum(logs) - minimum(logm)
    line1 = lines!(
        ax, [xmin, xmax], [ymin-a1, ymin + (xmax-xmin)],
        linewidth=.6, linestyle=:dot, color=:black
    )
    #~ Determine and plot line with b=2
    a2 = maximum(logs) - 2*maximum(logm)
    line2 = lines!(
        ax, [(ymin-a2)/2, (ymax-a2)/2], [ymin, ymax],
        linewidth=.6, linestyle=:dash, color=:black
    )

    #/ Add clarifying labels
    #~ compute rotation in "screen space"
    sx = 1 / (xmax - xmin)
    sy = 1 / (ymax - ymin)
    angle = atan(sy, sx)
    text!(1.5, 1.5+a1, rotation=angle, text=L"b=1", align=(:left,:top), fontsize=10)
    sx = 1 / (xmax - xmin)
    sy = 2 / (ymax - ymin)
    angle = atan(sy, sx)
    text!(1.2, 1.5*2+a2, rotation=angle, text=L"b=2", align=(:left,:bottom), fontsize=10)

    vlines!(ax, [0.], color=:gray, linestyle=:dash, linewidth=.5)
    # text!(0., ymax, rotation=π/2, text=L"Np=1", align=(:right,:bottom), color=:gray, fontsize=10)
    #~ Rare / common text
    text!(0.05,0.98, space=:relative, fontsize=9, text=L"\textrm{rare}", align=(:left,:top))
    text!(0.98,0.05, space=:relative, fontsize=9, text=L"\textrm{common}", align=(:right,:bottom))
    
    #~ Save
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module AFDPlotter
#/ End module
