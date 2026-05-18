#= Module to plot Figure 2 of main paper =#
module Figure2

using Meris
using DataFrames, DataFramesMeta, StatsBase
using CairoMakie, MakiePublication, LaTeXStrings
using Colors
using JLD2

include("./../scripts/module-scripts/macropatterns/taylor.jl")
using .Taylor

include("./../plot/plot-taylor.jl")
using .TaylorPlotter

include("./../plot/colors/shadetester.jl")
using .Shades: shades

### DATA PREPARATION ###
function _load_downsampled_category(category::AbstractString; downsampled_dir=Meris.DATADIR * "downsampled/")
    path = joinpath(downsampled_dir, "$(category).jld2")
    isfile(path) || error(
        "Missing grouped downsampled dataset: $path. " *
        "Run julia --project=julia-code julia-code/scripts/cli-scripts/downsample-dataset-groups.jl first."
    )

    d = JLD2.load(path)
    haskey(d, "ds_df") || error("Downsampled file does not contain ds_df: $path")

    df = d["ds_df"]
    select!(df, :class, :sample_id, :component_id, :counts, :nreads)

    return df
end

function prepare(;
        categories=["linguistic", "microbial", "social", "biology"],
        downsampled_dir=Meris.DATADIR * "downsampled/",
        outdir=Meris.DATADIR * "macro/taylor/",
    )
    mkpath(outdir)

    for category in categories
        @info "Loading downsampled $category data..."
        df = _load_downsampled_category(category; downsampled_dir=downsampled_dir)

        @info "Working $category data..."
        Taylor.compute(df; save=true, filename=joinpath(outdir, "$(category).jld2"))

        df = nothing
        GC.gc()
    end
end

### MAKE FIGURE 2 ###
function plot(;
        ext="pdf",
        big_limits=(-9, 0, -18, 0),
        small_limits=reverse([[-8, 0, -14, -3], [-8, 0, -14, -1], [-11, 0, -19, -1], [-12, 0, -22, 0]]),
        font_scale=1.5,
        height_scale=0.55,
        panel_colgap=6,
        small_rowgap=2,
        small_colgap=3,
        bases = [
            colorant"#1f77b4",  # blue
            colorant"#ff7f0e",  # orange
            colorant"#9467bd",  # purple
            colorant"#2ca02c",  # green
            colorant"#d62728"   # red
        ]
    )
    palette1 = shades(bases[1], 10)
    palette2 = shades(bases[2], 10)
    palette3 = shades(bases[3], 8)
    palette4 = vcat(shades(bases[4], 10)[1:7], shades(bases[5], 8))

    palettes = [palette1, palette2, palette3, palette4]

    fig = Figure(
        size = (TaylorPlotter.NATURE_DOUBLE_WIDTH_PT, height_scale * TaylorPlotter.NATURE_MAX_HEIGHT_PT),
        figure_padding = (4, 4, 4, 4)
    )
    TaylorPlotter.plot!(fig[1,1]; palettes=palettes, panel_start=1, font_scale=font_scale,
        big_limits=big_limits,
        small_limits=small_limits,
        panel_colgap=panel_colgap,
        small_rowgap=small_rowgap,
        small_colgap=small_colgap,
        center_data=false
        )
    save(Meris.FIGDIR * "fig2.$ext", fig, pt_per_unit=1)
    return fig
end

end
