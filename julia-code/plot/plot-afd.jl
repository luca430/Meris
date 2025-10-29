#= Module to plot Heap's law for some processes and/or data =#
#/ Start module
module AFDPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames
using StatsBase
using Distributions
using FHist

#/ Modules
import Meris.DATADIR as DATADIR
import Meris.MLEstimator as MLE
import Meris.StraightLine as SL

#################
### FUNCTIONS ###
function plot_afds(;
    DIRECTORIES = [ "rfc/", "lego/","bci.tree/"],
    LABELS = [L"\textrm{RFC}", L"\textrm{LEGO}", L"\textrm{BCI}"],
    nbins::Int = 27,
    DIR = DATADIR * "macro/afd/",
    BASEFILENAME = "z-values.csv",
    BASETLFILENAME = "tl-stats.csv",
    savefig = false,
    figname = true
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(1.5*width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"z", xlabelsize=12,
        ylabel=L"p(z)", ylabelsize=12,
        yscale=log10,
        limits=(-24,6,1e-8,1e1)
        # limits = (-5.,3.,0,0.6)
    )
    axtl = Axis(
        fig[1,2],
        xlabel=L"\textrm{mean}\;\mu^2", xlabelsize=12,
        ylabel=L"\textrm{variance}\;\sigma^2", ylabelsize=12,
        # limits=(-25,0,-25,0)
    )

    #/ Plot
    xtl = -25:0.1:0.0
    ytl = copy(xtl)
    lines!(axtl, xtl, ytl, linewidth=.8, color=:black, linestyle=(:dash,:dense))
    
    for i in eachindex(DIRECTORIES)
        zdf = CSV.read(DIR*DIRECTORIES[i]*BASEFILENAME, DataFrame)
        z = zdf[!,:z]
        zmin, zmax = extrema(zdf[!,:z])
        binedges = range(zmin, zmax, nbins)
        fh = FHist.Hist1D(z; counttype=Int, binedges=binedges, overflow=true) |> normalize
        #~ Scatter
        s = scatter!(
            ax, bincenters(fh), fh.bincounts, label=LABELS[i], markersize=4.5, strokewidth=1.
        )
        #~ try something
        # f = MLE.Burr
        # θstar = MLE.fit(f, exp.(z), [1.0,1.0,1.0])
        # xplot = exp.(-6:0.01:6)
        # yplot = xplot .* f.(xplot, Ref([0.1,1.5,3.0]))
        # lines!(ax, log.(xplot), yplot, linewidth=1., label=L"\textrm{Burr}")

        tldf = CSV.read(DIR*DIRECTORIES[i]*BASETLFILENAME, DataFrame)
        #~ Fix a straight line using York's method
        xs = -25:1.0:0.0
        
        
        x = log.(tldf[!,:meanfrequency].^2)
        y = log.(tldf[!,:varfrequency])
        wx = tldf[!,:noccurences] .* tldf[!,:meanfrequency].^2 ./ tldf[!,:varfrequency] ./ 4
        wy = (tldf[!,:noccurences] .- 1) ./ 2
        ρ = tldf[!,:thirdmomentfrequency] ./ (tldf[!,:noccurences] .* tldf[!,:meanfrequency] .* tldf[!,:varfrequency])
        ρ = nothing
        straightlinefit = SL.weightedyorkfit(x, y, wx, wy, ρ=ρ)
        # @info "fit" straightlinefit
        l = lines!(
            axtl, xs, straightlinefit.a .+ straightlinefit.b .* xs,
            linestyle=(:dash,:dense), linewidth=.8, color=colors[i]
        )
        #~ Scatter w.r.t. to their weight
        minwidth = .33
        maxwidth = 1.0
        strokewidth = (maxwidth - minwidth) .* wx./sum(wx) .+ minwidth
        #~ rescale so that they have slope 1
        # yrescaled = y .- straightlinefit.a
        stl = scatter!(axtl, x, y, markersize=4, strokewidth=strokewidth)
        # _mean, _var = tldf[!,:meanfrequency].^2, tldf[!,:varfrequency]
        # stl = scatter!(axtl, log.(_mean), log.(_var), markersize=4, strokewidth=.7)
    end

    
    axislegend(
        ax,
        position=:lt, labelsize=9, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

"Plot AFD of component systems of a very specific 'type'"
function plot_typedafd(;
    DIRECTORY = "macro/afd/lego/",
    LABELS = [L"\textrm{Star Wars}"],
    nbins::Int = 27,
    DIR = DATADIR,
    BASEFILENAME = "themed-z-values.csv",
    # BASETLFILENAME = "themed-tl-stats.csv",
    savefig = false,
    figname = true
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
        xlabel=L"z", xlabelsize=12,
        ylabel=L"p(z)", ylabelsize=12,
        yscale=log10,
        limits=(-25.0,5.0,1e-6,1e0)
    )
    
    #/ Load data
    zdf = CSV.read(DIR*DIRECTORY*BASEFILENAME, DataFrame)
    z = zdf[!,:z]
    zmin, zmax = extrema(zdf[!,:z])
    binedges = range(zmin, zmax, nbins)
    fh = FHist.Hist1D(z; counttype=Int, binedges=binedges, overflow=true) |> normalize
    #~ Scatter
    s = scatter!(ax, bincenters(fh), fh.bincounts, label=LABELS[begin], markersize=5)
    #~ do something stupid
    gamma = Distributions.fit_mle(Gamma, exp.(z))
    zgamma = -25.0:0.01:3.0
    gammaplot = exp.(zgamma) .* Distributions.pdf.(gamma, exp.(zgamma))
    lines!(ax, zgamma, gammaplot, color=:black, linestyle=(:dash,:dense), linewidth=1.)

    axislegend(
        ax,
        position=:lt, labelsize=9, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end


"Plot AFD of the Brown corpus"
function plot_brownafd(;
    methods = ["noreplace", "replace", "multinomial", "mvhypgeom"],
    nbins::Int = 19,
    DIR = DATADIR * "macro/afd/",
    FILENAME = "Nfixed_zvalues.csv",
    TLFILENAME = "Nfixed_stats.csv",
    savefig = false,
    figname = true
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(1.5*width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"z", xlabelsize=12,
        ylabel=L"p(z)", ylabelsize=12,
        # yscale=log10,
        limits = (-5.,3.,0,0.6)
    )
    axtl = Axis(
        fig[1,2],
        xlabel=L"\textrm{mean}\;\mu", xlabelsize=12,
        ylabel=L"\textrm{variance}\;\sigma^2", ylabelsize=12,
    )

    #/ Plot
    β = zeros(length(methods))
    a = zeros(length(methods))
    b = zeros(length(methods))
    μplot = -10.0:0.01:0.0
    for i in eachindex(methods)
        zdf = CSV.read(DIR*methods[i]*"_"*FILENAME, DataFrame)
        z = zdf[!,:z]
        zmin, zmax = extrema(zdf[!,:z])
        binedges = range(zmin, zmax, nbins)
        fh = FHist.Hist1D(z; counttype=Int, binedges=binedges, overflow=true) |> normalize
        #~ Compute β using Taylor's law
        tldf = CSV.read(DIR*methods[i]*"_"*TLFILENAME, DataFrame)
        c = mean(tldf[!,:meanfrequency] .^2 ./ tldf[!,:varfrequency])
        a[i], b[i] = CurveFit.linear_fit(
            log.(tldf[!,:meanfrequency]), log.(tldf[!,:varfrequency])
        )
        # ltl = lines!(axtl, μplot, af .+ bf .* μplot)
        β[i] = 2*c/(c-1)
        # stl = scatter!(axtl, 2.0 .* log.(tldf[!,:meanfrequency]), log.(tldf[!,:varfrequency]))
        stl = scatter!(axtl, log.(tldf[!,:meanfrequency]), log.(tldf[!,:varfrequency]))

        #~ Try to plot log-Lomax
        #~ Find the `b` param of the log-lomax by fitting [@TODO Use Taylor's law instead]
        # neglogll(p) = -sum(loglomax(z, p))
        # res = Optim.optimize(neglogll, 1e-6, 1e4)
        # btemp = Optim.minimizer(res)
        # b += btemp
        
        #~ Scatter
        s = scatter!(ax, bincenters(fh), fh.bincounts, label=L"\textrm{%$(methods[i])}")
    end

    βplot = mean(β)
    zplot = -5:0.01:3
    pplot = loglomax.(zplot, βplot)
    lines!(ax, zplot, exp.(pplot), linewidth=.8, color=:black, label=L"\textrm{Lomax}")

    #~ Plot a straight line with exponent β
    σplot(b) = b .* μplot
    lines!(axtl, μplot, mean(a) .+ σplot(mean(b)), linewidth=.8, color=:black)
    lines!(axtl, [-6.,-3.], [-14.,-11.], linewidth=.8, color=:black, linestyle=(:dash,:dense))
    text!(
        axtl, -5.25, -11.5; text=L"\sigma^2 \propto \mu",
        align=(:center,:center), rotation=π/4.5
    )

    
    axislegend(
        ax,
        position=:lt, labelsize=9, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

########################
### HELPER FUNCTIONS ###

#/ A ZOO OF DISTRIBUTIONS
function loggamma(z, α)
    return α*sqrt(trigamma(α)) .* z .+ α*digamma(α) .- exp.(z .* sqrt(trigamma(α)) .+ digamma(α)) .+ 0.5*log(trigamma(α)) .- loggamma(α)
end

"Logarithmic Lomax distribution"
function loglomax(z, b)
    s = sqrt(trigamma(1) + trigamma(b))
    m = digamma(1) - digamma(b)
    return log(s * b) .+ z .* s .+ m .- (b + 1) .* log.(1 .+ exp.(z .* s .+ m))
end

function lrln(z, σ)
    return -z .^ 2 ./ 2 .- log(sqrt(σ^2 * 2 * π))
end
##############################

end # module AFDPlotter
#/ End module
