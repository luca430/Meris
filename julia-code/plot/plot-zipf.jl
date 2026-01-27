#= Module to plot Zipf's law for some processes and/or data =#
#/ Start module
module ZipfPlotter

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
function plot_zipf(;
    BROWNDIR = DATADIR * "zipf/brown/",
    WIKIDIR = DATADIR * "zipf/wikitext/",
    LEGODIR = DATADIR * "zipf/lego/",
    TREEDIR = DATADIR * "zipf/bci.tree/",
    RFCDIR  = DATADIR * "zipf/rfc/",
    OTUDIR  = DATADIR * "zipf/otu/",
    reduce  = true,     #~ Reduces number of points plotted
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
        xlabel=L"\textrm{rank}\;\log_{10}r", xlabelsize=11,
        ylabel=L"\textrm{frequency}\;\log_{10} \nu", ylabelsize=11,
        limits=(0,6,-9,0)
    )
    #~ Define some functions for fitting
    logfun(x,p) = @. p[1] + p[2]*log(x)
    powfun(x,p) = @. p[1] + p[2]*x^p[3]
    
    #/ Load and plot wikitext
    wikitextfilename = WIKIDIR*"wikitext-103-zipf.jld2"
    _plot_zipf(
        ax, wikitextfilename;
        npoints = 2^5 - 1,
        marker=:dtriangle, strokecolor=:black,
        label=L"\textrm{Wikitext-103}"
    )

    #/ Load and plot Brown corpus
    brownfilename = BROWNDIR*"brown-zipf.jld2"
    _plot_zipf(
        ax, brownfilename;
        npoints = 2^5 - 1,
        marker=:diamond, strokecolor=:rebeccapurple,
        label=L"\textrm{Brown corpus}"
    )
    
    #/ Load and plot LEGO
    # legofilename = LEGODIR*"lego-zipf.jld2"
    # _plot_zipf(ax, legofilename; marker=:rect, strokecolor=colors[1], label=L"\textrm{LEGO}")
    # #/ Load and plot BCI trees
    # treefilename = TREEDIR*"bcitree-zipf.jld2"
    # _plot_zipf(ax, treefilename; marker=:xcross, strokecolor=colors[2], label=L"\textrm{BCI}")
    # #/ Load and plot OTUs
    # environments = ["gut1", "gut2", "seawater"]
    # for (i, env) in enumerate(environments)
    #     otufilename = OTUDIR*"otu-$(env)-zipf.jld2"
    #     _plot_zipf(
    #         ax, otufilename; marker=:circle, strokecolor=colors[2+i], label=L"\textrm{%$(env)}"
    #     )
    # end
    # #/ Load and plot RFC
    # rfcfilename = RFCDIR*"rfc-zipf.jld2"
    # _plot_zipf(
    #     ax, rfcfilename; marker=:utriangle, strokecolor=colors[3+length(environments)],
    #     label=L"\textrm{RFC}"
    # )
    
    #/ Add legend    
    axislegend(
        ax,
        position=:lb, labelsize=9, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

function plot_syntheticzipf(;
    ZIPFDIR = DATADIR * "zipf/synthetic/",
    filename = "synthetic-zipf-1.jld2",
    n = 31,
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
        fig[1,1], aspect=1,
        xlabel=L"\textrm{rank}\;\log_{10}r", xlabelsize=11,
        ylabel=L"\textrm{frequency}\;\log_{10}\nu", ylabelsize=11,        
        xminorticks=IntervalsBetween(4),
        yminorticks=IntervalsBetween(4),
        limits=(0,4,-6,0),
        xticks = range(0, 4, 5), yticks=range(-6,0,4)
    )

    #~ Load data
    db = JLD2.load(ZIPFDIR*filename)
    logr = log10.(db["ranks"])
    logf = log10.(db["freqs"])
    #~ subsample, as there's too many points otherwise
    dn = floor(Int, length(logr) / n)
    rplot = logr[1:dn:end]
    fplot = logf[1:dn:end]
    
    params = db["params"]
    #/ Scatter synthetic data
    scatter!(
        ax, rplot, fplot,
        color=(colors[3],0.7), strokecolor=:black, markersize=4, strokewidth=.4,
    )
    #~ Straight line, power law
    xmin, xmax, ymin, ymax = ax.limits[]
    xstart = 1.5
    ystart = -1.8
    yend = -5
    ζ = 1/params.γ
    xend = xstart + (ystart-yend)/ζ
    lines!(
        ax, [xstart, xend], [ystart, yend],
        color=:black, linestyle=(:dash,:dense), linewidth=.8
    )    
    text!(
        xstart+(xend-xstart)/2, ystart+(yend-ystart)/2, text=L"\propto r^{-\zeta(\gamma)}",
        rotation=_get_angle(ax, -ζ), align=(:center,:bottom), fontsize=12
    )

    #~ Save
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

########################
### HELPER FUNCTIONS ###
function _plot_zipf(
    ax::Axis,
    filename::String;
    npoints = 31,
    marker=:rect,
    strokecolor=:black,
    label=L"\textrm{LABEL}"
)
    #~ Extract relevant data, reduce points to `npoints`, and plot
    df = JLD2.load(filename)
    rank, freq = df["rank"], df["frequency"]

    if !isnothing(npoints)
        logmax = log10(rank[end])
        idxs = round.(Int, exp10.(range(0, logmax, length=npoints+1)))
        uidxs = unique(idxs)
        permutation = sortperm(rank)
        __rank = rank[permutation]
        __freq = freq[permutation]
        rank = __rank[uidxs]
        freq = __freq[uidxs]
    end
    
    scatter!(
        ax, log10.(rank), log10.(freq),
        marker=marker, markersize=4,
        color=:white, strokecolor=strokecolor, strokewidth=.42,
        label=label
    )
    nothing
end

end # module HeapPlotter
#/ End module
