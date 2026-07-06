# Reproduce figures from the checked-in analysis intermediates.

using Dates

root_dir = normpath(joinpath(@__DIR__, ".."))
julia_dir = joinpath(root_dir, "julia-code")
results_dir = joinpath(root_dir, "results")
figure_set = "main"

function parse_args!(args)
    global results_dir, figure_set

    for arg in args
        if startswith(arg, "--results-dir=")
            results_dir = abspath(last(split(arg, "="; limit=2)))
        elseif startswith(arg, "--figure-set=")
            figure_set = lowercase(last(split(arg, "="; limit=2)))
        elseif arg == "--help" || arg == "-h"
            println("""
            Usage:
              julia --project=julia-code codeocean/reproduce_figures.jl [options]

            Options:
              --figure-set=main|taylor|all   Figure group to render. Default: main
              --results-dir=DIR              Directory receiving copied outputs.
            """)
            exit(0)
        else
            error("Unknown argument: $arg")
        end
    end
end

parse_args!(ARGS)

fig_results_dir = joinpath(results_dir, "figures")
mkpath(fig_results_dir)

using Meris

function copy_figure(filename)
    src = joinpath(Meris.FIGDIR, filename)
    isfile(src) || error("Expected figure was not created: $src")
    cp(src, joinpath(fig_results_dir, filename); force=true)
    return nothing
end

function include_plot(filename)
    include(joinpath(julia_dir, "plot", filename))
end

function render_main_figures()
    include_plot("plot-fig2-A.jl")
    if !isfile(joinpath(julia_dir, "plot", "plot-tl-prediction-bins.jl"))
        @info "Using Figure2A as the shared TL-prediction bin helper expected by Figure2B"
        @eval const TLPredictionBinPlotter = Figure2A
    end
    include_plot("plot-fig2-B.jl")
    include_plot("plot-fig3.jl")

    @info "Rendering Figure 2A"
    Base.invokelatest(Figure2A.plot)
    copy_figure("fig2_A.pdf")

    @info "Rendering Figure 2B"
    Base.invokelatest(Figure2B.plot)
    copy_figure("fig2_B.pdf")

    @info "Rendering Figure 3"
    Base.invokelatest(Figure3.plot)
    copy_figure("fig3.pdf")
end

function render_taylor_figures()
    include_plot("plot-taylor.jl")
    include_plot("plot-tl-prediction.jl")

    @info "Rendering Taylor-law supplementary figures"
    Base.invokelatest(TaylorPlotter.plot_taylor_downsampled)
    copy_figure("taylor.pdf")

    Base.invokelatest(TaylorPlotter.plot_taylor_unbinned_downsampled)
    copy_figure("taylor-unbinned.pdf")

    Base.invokelatest(TaylorPlotter.plot_taylor_reference_downsampled)
    copy_figure("taylor-reference.pdf")

    Base.invokelatest(TaylorPlotter.plot_taylor_omega)
    copy_figure("taylor-omega.pdf")

    Base.invokelatest(TaylorPlotter.plot_taylor_omega_binned)
    copy_figure("taylor-omega-binned.pdf")

    Base.invokelatest(TaylorPlotter.plot_gutenberg_en_count_distribution)
    copy_figure("gutenberg-en-count-distribution.pdf")

    Base.invokelatest(TaylorPlotter.plot_single_domain_count_distributions)
    copy_figure("single-domain-count-distributions.pdf")

    Base.invokelatest(TaylorPlotter.plot_dataset_key)
    copy_figure("dataset-key.pdf")

    Base.invokelatest(TLPredictionPlotter.plot)
    copy_figure("tl-prediction-all.pdf")
end

if figure_set == "main"
    render_main_figures()
elseif figure_set == "taylor"
    render_taylor_figures()
elseif figure_set == "all"
    render_main_figures()
    render_taylor_figures()
else
    error("Unknown figure set: $figure_set. Use main, taylor, or all.")
end

open(joinpath(results_dir, "reproduction-summary.txt"), "w") do io
    println(io, "Reproduction completed at ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))
    println(io, "Figure set: ", figure_set)
    println(io, "Figures written to: ", fig_results_dir)
    println(io, "AIC tables written to: ", joinpath(results_dir, "tables"))
end
