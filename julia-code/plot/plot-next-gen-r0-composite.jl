#= Compose the beta scan and grid heatmap next-generation R0 figures. =#
module NextGenR0CompositePlot

using CairoMakie
using MakiePublication

using Meris

include("plot-next-gen-r0-beta-scan.jl")
using .NextGenR0BetaScanPlot

include("plot-next-gen-r0-grid.jl")
using .NextGenR0GridPlot

include("plot-next-gen-r0.jl")
using .NextGenR0Plot

const DEFAULT_BETA_INPUT = NextGenR0BetaScanPlot.DEFAULT_INPUT
const DEFAULT_GRID_INPUT = NextGenR0GridPlot.DEFAULT_INPUT
const DEFAULT_DISTRIBUTION_INPUT = NextGenR0Plot.DEFAULT_INPUT
const DEFAULT_OUTDIR = joinpath(Meris.FIGDIR, "next-gen")
function parse_args(args)
    options = Dict(
        "beta-input" => DEFAULT_BETA_INPUT,
        "grid-input" => DEFAULT_GRID_INPUT,
        "distribution-input" => DEFAULT_DISTRIBUTION_INPUT,
        "outdir" => DEFAULT_OUTDIR,
        "basename" => "otu-gut1-r0-composite",
        "probability-gamma" => "0.35",
        "distribution-nbins" => "28",
    )

    for arg in args
        if arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code julia-code/plot/plot-next-gen-r0-composite.jl [options]

            Options:
              --beta-input=PATH        Input JLD2 beta-scan file. Default: $(DEFAULT_BETA_INPUT)
              --grid-input=PATH        Input JLD2 grid file. Default: $(DEFAULT_GRID_INPUT)
              --distribution-input=PATH
                                       Input JLD2 R0 distribution file. Default: $(DEFAULT_DISTRIBUTION_INPUT)
              --outdir=PATH            Output directory. Default: $(DEFAULT_OUTDIR)
              --basename=NAME          Output filename stem. Default: otu-gut1-r0-composite
              --probability-gamma=G    Heatmap display gamma. Default: 0.35
              --distribution-nbins=N   Number of bins for R0 distribution. Default: 28
            """)
            exit(0)
        elseif startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "="; limit=2)
            haskey(options, key) || error("Unknown option: --$key")
            options[key] = value
        else
            error("Unknown argument: $arg")
        end
    end

    return (
        beta_input = options["beta-input"],
        grid_input = options["grid-input"],
        distribution_input = options["distribution-input"],
        outdir = options["outdir"],
        basename = options["basename"],
        probability_gamma = parse(Float64, options["probability-gamma"]),
        distribution_nbins = parse(Int, replace(options["distribution-nbins"], "_" => "")),
    )
end

function plot(;
    beta_input::AbstractString=DEFAULT_BETA_INPUT,
    grid_input::AbstractString=DEFAULT_GRID_INPUT,
    distribution_input::AbstractString=DEFAULT_DISTRIBUTION_INPUT,
    outdir::AbstractString=DEFAULT_OUTDIR,
    basename::AbstractString="otu-gut1-r0-composite",
    probability_gamma::Real=0.35,
    distribution_nbins::Int=28,
    savefig::Bool=true,
)
    set_theme!(MakiePublication.theme_acs(; ishollowmarkers=[true, true]))

    fig = Figure(; size=(590, 520), figure_padding=(70, 34, 24, 30))

    Label(
        fig[0, 1],
        "a.";
        fontsize=23,
        font=:bold,
        halign=:left,
        valign=:bottom,
        padding=(0, 0, 0, 0),
    )
    Label(
        fig[0, 2],
        "b.";
        fontsize=23,
        font=:bold,
        halign=:left,
        valign=:bottom,
        padding=(0, 0, 0, 0),
    )
    left = GridLayout(fig[1, 1])
    right = GridLayout(fig[1, 2])
    bottom_wrapper = GridLayout(fig[3, 1:2])
    bottom = GridLayout(bottom_wrapper[1, 1])

    Label(
        bottom_wrapper[0, 1],
        "c.";
        fontsize=23,
        font=:bold,
        halign=:left,
        valign=:bottom,
        padding=(0, 0, 0, 0),
    )

    beta_panel = NextGenR0BetaScanPlot.plot!(
        left;
        input=beta_input,
        labelsize=20,
        ticklabelsize=13,
        textsize=15,
        markersize=7,
        linewidth_scale=1.05,
        legend_position=(0.998, 0.002),
        legend_margin=(0, 1, 0, 0),
        legend_patchlabelgap=1,
        axis_kwargs=(;
            spinewidth=1.1,
            xtickalign=0,
            ytickalign=0,
            xminortickalign=0,
            yminortickalign=0,
            xticksize=4,
            yticksize=4,
            xminorticksize=3,
            yminorticksize=3,
            xtickwidth=1.2,
            ytickwidth=1.2,
            xminortickwidth=0.9,
            yminortickwidth=0.9,
            xticklabelfont="TeX Gyre Heros",
            yticklabelfont="TeX Gyre Heros",
        ),
    )

    grid_panel = NextGenR0GridPlot.plot!(
        right;
        input=grid_input,
        probability_gamma=probability_gamma,
        labelsize=20,
        ticklabelsize=13,
        textsize=15,
        linewidth_scale=1.05,
        colorbar_width=8,
        colorbar_gap=0,
        colorbar_spinewidth=1.1,
        axis_kwargs=(;
            aspect=AxisAspect(1),
            spinewidth=1.1,
            xtickalign=0,
            ytickalign=0,
            xminortickalign=0,
            yminortickalign=0,
            xticksize=4,
            yticksize=4,
            xminorticksize=3,
            yminorticksize=3,
            xtickwidth=1.2,
            ytickwidth=1.2,
            xminortickwidth=0.9,
            yminortickwidth=0.9,
            xticklabelfont="TeX Gyre Heros",
            yticklabelfont="TeX Gyre Heros",
        ),
    )

    xlims!(beta_panel.axis, 5e-6, 1e-1)
    ylims!(beta_panel.axis, 0.0, 1.0)
    beta_panel.axis.yticks = ([0.0, 0.5, 1.0], ["0.0", "0.5", "1.0"])

    xlims!(grid_panel.axis, 0.0, 3.0)
    ylims!(grid_panel.axis, -5.0, -1.0)
    grid_panel.axis.xticks = ([0.0, 1.0, 2.0, 3.0], ["0", "1", "2", "3"])
    grid_panel.axis.yticks = ([-5.0, -4.0, -3.0, -2.0, -1.0], ["−5", "−4", "−3", "−2", "−1"])

    distribution_panel = NextGenR0Plot.plot!(
        bottom;
        input=distribution_input,
        nbins=distribution_nbins,
        labelsize=20,
        ticklabelsize=13,
        linewidth_scale=1.25,
        legend_labelsize=14,
        axis_kwargs=(;
            spinewidth=1.1,
            xtickalign=0,
            ytickalign=0,
            xminortickalign=0,
            yminortickalign=0,
            xticksize=5,
            yticksize=5,
            xminorticksize=4,
            yminorticksize=4,
            xtickwidth=1.5,
            ytickwidth=1.5,
            xminortickwidth=1.1,
            yminortickwidth=1.1,
            xticklabelfont="TeX Gyre Heros",
            yticklabelfont="TeX Gyre Heros",
        ),
    )
    xlims!(distribution_panel.axis, -4.0, 3.0)
    ylims!(distribution_panel.axis, 0.0, 1.5)
    distribution_panel.axis.xticks = (
        collect(-4.0:1.0:3.0),
        ["−4", "−3", "−2", "−1", "0", "1", "2", "3"],
    )
    distribution_panel.axis.yticks = ([0.0, 0.5, 1.0, 1.5], ["0.0", "0.5", "1.0", "1.5"])

    colsize!(right, 1, Fixed(112))
    colsize!(right, 2, Fixed(18))
    colgap!(right, -2)

    colsize!(bottom_wrapper, 1, Fixed(430))
    rowsize!(bottom_wrapper, 0, Fixed(18))
    rowgap!(bottom_wrapper, 0)

    rowsize!(fig.layout, 0, Fixed(30))
    rowsize!(fig.layout, 1, Fixed(128))
    rowsize!(fig.layout, 2, Fixed(0))
    rowsize!(fig.layout, 3, Fixed(154))
    colsize!(fig.layout, 1, Fixed(260))
    colsize!(fig.layout, 2, Fixed(155))
    rowgap!(fig.layout, 5)
    colgap!(fig.layout, -12)

    if savefig
        mkpath(outdir)
        pdf_file = joinpath(outdir, "$basename.pdf")
        png_file = joinpath(outdir, "$basename.png")
        CairoMakie.save(pdf_file, fig, pt_per_unit=1)
        CairoMakie.save(png_file, fig, px_per_unit=3)
        @info "Saved R0 composite plot" pdf=pdf_file png=png_file
    end

    return fig
end

end # module NextGenR0CompositePlot

if abspath(PROGRAM_FILE) == @__FILE__
    options = NextGenR0CompositePlot.parse_args(ARGS)
    NextGenR0CompositePlot.plot(;
        beta_input=options.beta_input,
        grid_input=options.grid_input,
        distribution_input=options.distribution_input,
        outdir=options.outdir,
        basename=options.basename,
        probability_gamma=options.probability_gamma,
        distribution_nbins=options.distribution_nbins,
    )
end
