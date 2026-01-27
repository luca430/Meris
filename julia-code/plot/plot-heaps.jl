#= Module to plot Heaps' law for some processes and/or data =#
#/ Start module
module HeapsPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings
using LsqFit
using JLD2
using StatsBase

#/ Modules
import Meris.DATADIR as DATADIR     #~ Data with (parsed/analyzed) results, for plotting
import Meris.ARXIVDIR as RARXIVDIR  #~ Directory with raw arXiv data

#################
### FUNCTIONS ###
function plot_heaps(;
    RFCDIR = DATADIR * "heaps/rfc/",
    PYDIR = DATADIR * "heaps/pitmanyor/",
    # BOOKDIR = DATADIR * "heap/books/",
    LEGODIR = DATADIR * "heaps/lego/",
    TREEDIR = DATADIR * "heaps/bci.tree/",
    OTUDIR  = DATADIR * "heaps/otu/",
    # ARXIVDIR = DATADIR * "heap/arxiv/",
    # ARXIVCATEGORIES = ["q-bio.QM"], #readlines(RARXIVDIR*"categories.txt"),
    rescale = true,
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
        xlabel=L"\textrm{total observation length}\;\log_{10}N", xlabelsize=11,
        ylabel=L"\textrm{vocabulary size}\;\log_{10} V(N)", ylabelsize=11,
        limits=(1,6,0,5)
    )
    
    #~ Define some functions for fitting
    logfun(x,p) = @. p[1] + p[2]*log(x)
    powfun(x,p) = @. p[1] + p[2]*x^p[3]

    lines!(ax, [1,6], [1,6], color=:black, linestyle=(:dash,:dense), linewidth=.6)
    
    # SYNTHETIC DATA / SIMULATIONS
    # for i in eachindex(αv)
    #     αs = round(αv[i], digits=1)
    #     # if αv[i] == 0.0
    #     #     # ~ Fit a line
    #     #     Nmin = 1
    #     #     Nmax = 6
    #     #     Nfit = exp10.(range(Nmin, Nmax, 128))
    #     #     llog = lines!(
    #     #         ax, log10.(Nfit), log10.(log.(Nfit)), color=colors[1],
    #     #         linewidth=1., label=L"\alpha=%$(αs)\;\textrm{(DP)}"
    #     #     )
    #     filename = "pitmanyor-heaps_a$(αs).jld2"
    #     db = JLD2.load(PYDIR*filename)
    #     N, V = db["N"], db["V"]

    #     # ~ Fit a line
    #     Nmin = 1
    #     Nmax = 6
    #     Nfit = exp10.(range(Nmin, Nmax, 128))
    #     #~ Idenfity index in data that is closest to Nmin
    #     idx = argmin(abs.(log10.(N) .- Nmin))
    #     pyfit = LsqFit.curve_fit(powfun, N[idx:end], V[idx:end], [1.,1.,.5])
    #     if rescale
    #         V = @. (V - pyfit.param[1]) / pyfit.param[2]
    #     end
    #     l = lines!(
    #         ax, log10.(Nfit), log10.(Nfit.^αv[i]), color=colors[i+3],
    #         linewidth=1., label=L"\alpha=%$(αs)\;\textrm{(PY)}"
    #     )            
    #     spow = scatter!(
    #         ax, log10.(N[idx:end]), log10.(V[idx:end]),
    #         marker=:circle, markersize=4, color=:white, strokecolor=colors[i+3]
    #     )
    # end
    
    #/ DATA
    #/ Load vocabulary for books
    # bookfilename = "stone-vocabsize.jld2"
    # bookdb = JLD2.load(BOOKDIR*"chinese/"*bookfilename)
    # Vbook, Nbook = bookdb["V"], bookdb["N"]
    # Vbook = dropdims(mean(Vbook, dims=2), dims=2)
    # #~ Fit a log-curve between Nmin and Nmax
    # Nmin = 3.25
    # Nmax = log10.(Nbook[end]) + 1
    # Nfit = exp10.(range(Nmin, Nmax, 128))
    # #~ Idenfity index in data that is closest to Nmin
    # idx = argmin(abs.(log10.(Nbook) .- Nmin))
    
    # bookfit = LsqFit.curve_fit(logfun, Nbook[idx:end], Vbook[idx:end], [1.,1.])
    # if rescale
    #     Vbook = @. (Vbook - bookfit.param[1]) / bookfit.param[2]
    # end
    # #~ scatter plot raw data
    # bookscat = scatter!(
    #     ax, log10.(Nbook), log10.(Vbook), marker=:rect, markersize=4,
    #     label=L"\textrm{Chinese book}", color=:white, strokecolor=colors[1]
    # )    
    # #~ Some text
    # ct = text!(
    #     ax, 5.5, 0.8, color=:black,
    #     text = L"\propto \log N", fontsize=10, align=(:center,:top)
    # )

    Nmax = 7
    
    #/ RFC documents
    rfcfilename = RFCDIR*"rfc-heaps.jld2"
    rfcplot = _scatter(ax, rfcfilename, colors[1]; label=L"\textrm{RFC}", Nmax=Nmax)
    #/ LEGO sets
    legofilename = LEGODIR*"lego-heaps.jld2"
    legoplot = _scatter(ax, legofilename, colors[2]; label=L"\textrm{LEGO}", Nmax=Nmax)
    #/ BCI trees
    treefilename = TREEDIR*"bci.tree-heaps.jld2"
    treeplot = _scatter(ax, treefilename, colors[3]; label=L"\textrm{BCI}", Nmax=Nmax)
    #/ OTUs
    otufilename = OTUDIR*"otu-gut1-heaps.jld2"
    otuplot = _scatter(ax, otufilename, colors[4]; label=L"\textrm{OTU}", Nmax=Nmax)

    # #/ Load vocabulary for LEGO
    # legofilename = "lego-vocabsize.jld2"
    # legodb = JLD2.load(LEGODIR*legofilename)
    # Vlego, Nlego = legodb["V"], legodb["N"]
    # Vlego = dropdims(mean(Vlego, dims=2), dims=2)
    # #~ Idenfity index in data that is closest to Nmin
    # Nmin = 1.0
    # Nfit = exp10.(range(Nmin, Nmax, 128))
    # idx = argmin(abs.(log10.(Nlego) .- Nmin))
    # legofit = LsqFit.curve_fit(powfun, Nlego[idx:end], Vlego[idx:end], [1.,1.,.5])
    # if rescale
    #     Vlego = @. (Vlego - legofit.param[1]) / legofit.param[2]
    # end
    # #~ scatter plot raw data
    # l = lines!(
    #     ax, log10.(Nfit), log10.(Nfit.^legofit.param[3]),
    #     linewidth=1., color=colors[2], linestyle=(:dashdot,:dense)
    # )
    # legoscat = scatter!(
    #     ax, log10.(Nlego), log10.(Vlego), marker=:rect, markersize=4,
    #     label=L"\textrm{LEGO}", color=:white, strokecolor=colors[2]
    # )
    # #~ Some text
    # αlego = round(legofit.param[3], digits=2)
    # lt = text!(
    #     ax, 5.5, 2.5, color=:black, rotation=0.175,
    #     text = L"\propto N^{%$(αlego)}", fontsize=10, align=(:center,:top)
    # )

    # #/ Load vocabulary for OTUs
    # otufilename = "otu-gut1-vocabsize.jld2"
    # otudb = JLD2.load(OTUDIR*otufilename)
    # Votu, Notu = otudb["V"], otudb["N"]
    # Votu = dropdims(mean(Votu, dims=2), dims=2)
    # #~ Idenfity index in data that is closest to Nmin
    # Nmin = 1.0
    # Nfit = exp10.(range(Nmin, Nmax, 128))
    # idx = argmin(abs.(log10.(Notu) .- Nmin))
    # otufit = LsqFit.curve_fit(powfun, Notu[idx:end], Votu[idx:end], [1.,1.,.5])
    # if rescale
    #     Votu = @. (Votu - otufit.param[1]) / otufit.param[2]
    # end
    # #~ scatter plot raw data
    # l = lines!(
    #     ax, log10.(Nfit), log10.(Nfit.^otufit.param[3]),
    #     linewidth=1., color=colors[3], linestyle=(:dashdot,:dense)
    # )
    # otuscat = scatter!(
    #     ax, log10.(Notu), log10.(Votu), marker=:rect, markersize=4,
    #     label=L"\textrm{OTU-gut1}", color=:white, strokecolor=colors[3]
    # )
    # #~ Some text
    # αotu = round(otufit.param[3], digits=2)
    # ot = text!(
    #     ax, 5.5, 3.5, color=:black, rotation=0.275,
    #     text = L"\propto N^{%$(αotu)}", fontsize=10, align=(:center,:bottom)
    # )

    # #/ Load vocabulary for arXiv papers [of a specific category]
    # for CATEGORY in ARXIVCATEGORIES
    #     arxivfilename = "arxiv-$(CATEGORY)-vocabsize.jld2"
    #     arxivdb = JLD2.load(ARXIVDIR*arxivfilename)
    #     arxivdf, meanarxivdf = arxivdb["raw"], arxivdb["average"]
    #     Narxiv, mNarxiv = arxivdf[!,:documentsize], meanarxivdf[!,:documentsize]
    #     Varxiv, mVarxiv = arxivdf[!,:vocabularysize], meanarxivdf[!,:vocabularysize]
    #     idx = argmin(abs.(log10.(Narxiv) .- Nmin))
    #     arxivfit = LsqFit.curve_fit(powfun, Narxiv[idx:end], Varxiv[idx:end], [1.,1.,.5])
    #     if rescale
    #         Varxiv = @. (Varxiv - arxivfit.param[1]) / arxivfit.param[2]
    #         mVarxiv = @. (mVarxiv - arxivfit.param[1]) / arxivfit.param[2]
    #     end
    #     #~ line plot
    #     l = lines!(
    #         ax, log10.(Nfit), log10.(Nfit.^arxivfit.param[3]),
    #         linewidth=1., color=colors[4], linestyle=(:dashdot,:dense)
    #     )
    #     #~ scatter plot raw data
    #     categorylabel = replace(CATEGORY, "-" => "‐")
    #     arxivscat = scatter!(
    #         ax, log10.(mNarxiv), log10.(mVarxiv), marker=:rect, markersize=4,
    #         # label=L"\textrm{arXiv-q‐bio.PE}",
    #         label=L"\textrm{arXiv-%$(categorylabel)}",
    #         color=:white, strokecolor=colors[4]
    #     )
    # end
    
    #/ Add legend    
    axislegend(
        ax,
        position=:lt, labelsize=10, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

function plot_syntheticheaps(;
    HEAPSDIR = DATADIR * "heaps/synthetic/",
    filename = "synthetic-heaps.jld2",
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
    
    #/ Plot Zipf's law for synthetic data
    ax = Axis(
        fig[1,1],
        xlabel=L"\textrm{sample size}\;\log_{10}N", xlabelsize=11,
        ylabel=L"\textrm{vocab. size}\;\log_{10}V", ylabelsize=11,        
        xminorticks=IntervalsBetween(4),
        yminorticks=IntervalsBetween(4),
        limits=(1,11,0,7),
        xticks=[1,3,5,7,9,11], yticks=[0,2,4,6]
    )

    #~ Load data
    db = JLD2.load(HEAPSDIR*filename)
    params = db["params"]
    logN = log10.(db["N"])
    logV = log10.(db["V"] ./ params.K)
    #/ Scatter synthetic data
    scatter!(
        ax, logN, logV,
        color=(colors[2],0.7), strokecolor=:black, markersize=4, strokewidth=.4
    )
    #~ Straight line, power law
    xmin, xmax, ymin, ymax = ax.limits[]
    # xstart = 1.0
    # xend = 2.5
    # ystart = 0.5
    # lines!(
    #     ax, [xstart, xend], [ystart, ystart + (xend-xstart)],
    #     linestyle=(:dash,:dense), color=:black, linewidth=.8
    # )
    η = min(params.γ, 1)

    xstart = 3.25
    xend = 9.5
    ystart = 2.05
    lines!(
        ax, [xstart,xend], [ystart, ystart + (xend-xstart)*η],
        # [ymin+2.5,(6.5-3)*β-ymin+2.5],
        linestyle=(:dash,:dense), color=:black, linewidth=.8
    )
    hlines!(ax, [log10.(params.S)], linestyle=:dot, color=:gray, linewidth=.5)    
    #~ Slope 1
    # xs = 1.8    # x-position, chosen by eye
    # ys = 1.2    # y-position, chosen by eye
    # text!(xs, ys, rotation=_get_angle(ax, 1), text=L"\propto\!N", align=(:center,:top), fontsize=10)
    #~ Slope η
    xs = 6.6
    ys = 3.4
    text!(
        xs, ys, rotation=_get_angle(ax, η), text=L"\propto\!N^{\eta(\gamma)}",
        align=(:center,:top), fontsize=12
    )

    #~ Convergence to system size S
    text!(
        1.2, log10(params.S), text=L"V(N)\rightarrow S", fontsize=10, color=:grey,
        align=(:left,:bottom)
    )

    #~ Save
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig    
end

########################
### HELPER FUNCTIONS ###
function _scatter(ax, FILENAME, strokecolor;
    Nmin=nothing, Nmax=7,
    rescale::Bool=false,
    label=nothing
)
    db = JLD2.load(FILENAME)
    N, V = db["N"], db["V"]
    # N = N / N[end]
    #/ Fit a straight line
    fitres = fit_line(N, V; Nmin=Nmin, Nmax=Nmax)
    #~ rescale the raw data when desired
    if rescale
        V = @. (V - fitres.params[1]) / fitres.params[2]
    end
    #~ plot raw data
    scat = scatter!(
        ax, log10.(N), log10.(V), markersize=4,
        color=:white, strokecolor=(strokecolor, 0.6), label=label
    )
    return (; s=scat)
    #~ overlay fitted line
    # line = lines!(
    #     ax, fitres.N, fitres.V,
    #     linewidth=1., color=strokecolor, linestyle=(:dashdot,:dense)
    # )
    # return (; s=scat, l=line)
end

function _get_angle(ax, γ)
    (xmin, xmax, ymin, ymax) = ax.limits[]
    sx = 1 / (xmax - xmin)
    sy = γ / (ymax - ymin)
    angle = atan(sy, sx)
    return angle
end

function fit_line(N, V; Nmin=nothing, Nmax=nothing)
    #~ Idenfity index in data that is closest to Nmin
    Nmin = 1
    isnothing(Nmax) && (Nmax = 1+log10(maximum(N)))
    Nfit = exp10.(range(Nmin, Nmax, 128))
    idx = argmin(abs.(log10.(N) .- Nmin))

    fit = LsqFit.curve_fit(powfun, N[idx:end], V[idx:end], [1.,1.,.5])
    return (N=log10.(Nfit), V=log10.(Nfit.^fit.param[3]), params=fit.param)
end

powfun(x,p) = @. p[1] + p[2]*x^p[3]

end # module HeapPlotter
#/ End module
