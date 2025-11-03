#= Simple module to check pure Pareto tails of a Burr distribution =#
#/ Start module
module BurrPlotter

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
    DATADIR = Meris.DATADIR * "heavy-tails/",
    FILENAME = "burr-power-law.jld2",
    plottail = false
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    #/ Load data
    d = JLD2.load(DATADIR*FILENAME)
    params = d["params"]
    fh = d["fh"]
    paretofit = d["paretofit"]
    tpdf = d["p"]

    logxmin, logxmax = log10(tpdf.x[begin]), log10(tpdf.x[end])

    #/ Make figure
    width = .9 * 246
    height = 3*width / 4.37
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"\log_{10}(x)", xlabelsize=11,
        ylabel=L"\textrm{rescaled}\;\log_{10}(p)", ylabelsize=11,
        limits=(logxmin,logxmax,-6,1)
    )

    #/ Plot    
    lines!(
        ax, log10.(tpdf.x), log10.(log(10) .* tpdf.x .* tpdf.y),
        linewidth=1., label=L"\textrm{Burr}"
    )
    scatter!(ax, fh.x, log10.(fh.y), label=L"\textrm{histogram}")

    #~ Plot power law tail
    if plottail
        fhtail = d["fhtail"]
        scatter!(
            ax, fhtail.x, log10.(fhtail.y), label=L"\textrm{tail histogram}",
            strokecolor=:grey
        )
    end

    #~ Plot pure Pareto fitted on the tail
    #  note: γ ← γ - 1, due to different conventions
    xpareto = collect(range(log10.(paretofit.xmin), logxmax, 32))
    γ = paretofit.γ - 1
    FPareto = Meris.ParetoLike.Paretocdf(γ; xmin=paretofit.xmin)
    ZPareto = FPareto(exp10(logxmax)) - FPareto(paretofit.xmin)
    #~ Scale it appropriately
    scale = Meris.ParetoLike.Burrpdf(paretofit.xmin, params.c, params.α, params.λ)
    scale = scale / Meris.ParetoLike.Paretopdf(paretofit.xmin, γ; xmin=paretofit.xmin)
    ypareto = Meris.ParetoLike.Paretopdf(exp10.(xpareto), γ; xmin=paretofit.xmin) ./ ZPareto
    ypareto = ypareto .* scale
    ypareto = log10.(log(10) .* exp10.(xpareto) .* ypareto)
    lines!(
        ax, xpareto, ypareto,
        color=:black, linewidth=1., linestyle=(:dash,:dense),
        label=L"\textrm{Pareto}"
    )
    vlines!(
        ax, [log10(paretofit.xmin)], ymin=-6, ymax=1,
        color=:black, linewidth=.8, linestyle=(:dot,:dense)
    )
    # add some text
    γstr = round(γ, digits=2)
    xtext = xpareto[24]
    ytext = ypareto[searchsortedfirst(xpareto, xtext)]
    tx = text!(
        xtext, ytext - 0.5,
        text=L"\propto x^{%$(γstr)}", fontsize=10,
        rotation = -π/4, align=(:center,:top)
    )

    #/ Add legend
    axislegend(ax, position=:lt, patchsize=(8,8), rowgap=0., labelsize=10, padding=0)
    
    return fig
end

end # module BurrPlotter
#/ End module
