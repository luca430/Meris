#= Module to plot Heap's law for some processes and/or data =#
#/ Start module
module AFDPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta
using StatsBase
using Distributions
using FHist

using JLD2

#/ Modules
import Meris.DATADIR as DATADIR
import Meris.StraightLine as SL

#################
### FUNCTIONS ###
function plot_afds(;
    DIRECTORIES = [ "rfc/", "lego/"], #, "bci.tree/"],
    LABELS = [L"\textrm{RFC}", L"\textrm{LEGO}", L"\textrm{BCI}"],
    nbins::Int = 27,
    DIR = DATADIR * "macro/afd/",
    BASEFILENAME = "z-values.csv",
    BASETLFILENAME = "tl-stats.csv",
    rescale = false,
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
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=12,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=12,
        limits=(-15,0,-25,0)
    )

    #/ Plot
    xtl = -20:0.1:0.0
    ytl = copy(xtl).*2
    lines!(axtl, xtl, ytl, linewidth=1, color=:black, linestyle=(:solid,:dense))
    
    for i in eachindex(DIRECTORIES)
        zdf = CSV.read(DIR*DIRECTORIES[i]*BASEFILENAME, DataFrame)

        #/ Abundance fluctuation distribution
        z = zdf[!,:z]
        zmin, zmax = extrema(zdf[!,:z])
        binedges = range(zmin, zmax, nbins)
        fh = FHist.Hist1D(z; counttype=Int, binedges=binedges, overflow=true) |> normalize
        #~ Scatter
        s = scatter!(
            ax, bincenters(fh), fh.bincounts, label=LABELS[i], markersize=4.5, strokewidth=1.
        )

        #/ Taylor's law
        tldf = CSV.read(DIR*DIRECTORIES[i]*BASETLFILENAME, DataFrame)
        m = tldf[!,:meanfrequency]
        s = tldf[!,:varfrequency]
        x = log.(m)
        y = log.(s)

        #~ Compute rescaled moments using the occupancy        
        tldf = @chain tldf begin 
            @transform(:omeanfrequency = :meanfrequency .* :occupancy)
            @transform(:ovarfrequency = :varfrequency .+ (1 .- :occupancy) .* :meanfrequency.^2)
            @transform(:ovarfrequency = :ovarfrequency .* :occupancy)
        end
        #~ Use them if `rescale=true`
        if rescale
            x = log.(tldf[!,:omeanfrequency])
            y = log.(tldf[!,:ovarfrequency])
        end
        
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
        straightlinefit = SL.weightedyorkfit(x, y, wx, wy, ρ=logρ)
        @info "fit" DIRECTORIES[i] straightlinefit
        #~ Reshuffle before plotting so it goes through (0,0)        
        xs = -20:1.0:0.0
        bs = round(straightlinefit.b, sigdigits=3)
        σs = round(straightlinefit.σb, sigdigits=2)
        blabel = L"b = %$(bs)\,(%$(σs))"
        l = lines!(
            axtl, xs, straightlinefit.b .* xs, label=blabel,
            linestyle=(:dash,:dense), linewidth=.8, color=colors[i]
        )
        #~ Scatter w.r.t. to their weight
        minwidth = .33
        maxwidth = 1.0
        strokewidth = (maxwidth - minwidth) .* wx./sum(wx) .+ minwidth
        #~ shift so that they go through the origin
        yshifted = y .- straightlinefit.a
        stl = scatter!(axtl, x, yshifted, markersize=4, strokewidth=strokewidth)
        # _mean, _var = tldf[!,:meanfrequency].^2, tldf[!,:varfrequency]
        # stl = scatter!(axtl, log.(_mean), log.(_var), markersize=4, strokewidth=.7)
    end

    
    axislegend(
        ax,
        position=:lt, labelsize=9, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    axislegend(
        axtl,
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

function plot_syntheticafd(;
    AFDDIR = DATADIR * "afd/synthetic/",
    filename = "synthetic-afd.jld2",
    nbins::Int = 31,
    savefig = false,
    figname = nothing
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
        xlabel=L"\log_{10}\,z", xlabelsize=11,
        ylabel=L"\textrm{pdf}\;\log_{10}\,p(z)", ylabelsize=11,
        xminorticks=IntervalsBetween(4),
        yminorticks=IntervalsBetween(4),
        limits=(-8,4,-3,0)
    )

    #~ Load data
    db = JLD2.load(AFDDIR*filename)
    z = db["z"]
    params = db["params"]

    #~ Fit histogram
    zmin, zmax = extrema(z)
    binedges = range(zmin, zmax, nbins)
    fh = FHist.Hist1D(z; binedges=binedges, overflow=true) |> normalize
    zs = bincenters(fh)
    ps = fh.bincounts

    #@TODO Compute a Gamma from system parameters instead
    xmin, xmax, ymin, ymax = ax.limits[]
    G = Distributions.fit_mle(Gamma, exp.(z))
    xgamma = exp.(range(xmin, xmax, 127))
    ygamma = xgamma .* Distributions.pdf.(G, xgamma)
    
    #/ Scatter synthetic data
    scatter!(
        ax, zs, log10.(ps),
        color=:white, strokecolor=:black, markersize=4, strokewidth=.4
    )

    lines!(
        ax, log.(xgamma), log10.(ygamma), label=L"\textrm{Gamma}",
        color=:black, linestyle=(:dash,:dense), linewidth=.8
    )
    #/ Legend and labels
    axislegend(
        ax,
        position=:cb, labelsize=8, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    text!(
        0.05,0.975, text=L"Np \gg 1", fontsize=10, align=(:left,:top), space=:relative,
        color=:grey
    )
    
    #~ Save
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
