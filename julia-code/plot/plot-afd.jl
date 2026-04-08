#= Module to plot AFD panels in a Figure 2-style 2x2 layout =#
module AFDPlotter

using Meris
using DataFrames
using CairoMakie, MakiePublication, LaTeXStrings
using Colors, ColorTypes
using FileIO, ImageTransformations

include("./../plot/colors/shadetester.jl")
using .Shades: shades

const ICONDIR = Meris.FIGDIR .* "icons"
const MM_TO_PT = 72.0 / 25.4
const NATURE_DOUBLE_WIDTH_PT = 183.0 * MM_TO_PT
const NATURE_MAX_HEIGHT_PT = 170.0 * MM_TO_PT
const NATURE_AXIS_LABEL_PT = 7
const NATURE_TICK_PT = 6

function ax_afd(ax, df, colors, markers; nbins::Int=25, occ::Float64=0.999, min_points::Int=50)
    classes = unique(df.class)
    xs = Float64[]
    ys = Float64[]

    for (i, class) in enumerate(classes)
        sdf = df[df.class .== class, :]
        afd = Meris.AFD.compute(sdf, :component_id; maxfrequency=1.0, minoccupancy=occ)

        length(afd.z) < min_points && continue

        x, y = Meris.DataTools.make_hist(afd.z; nbins=nbins)

        scatter!(
            ax, x, y;
            label=string(class),
            color=:white,
            strokecolor=colors[i],
            marker=markers[mod1(i, length(markers))],
            markersize=8,
            strokewidth=0.8,
        )

        append!(xs, x)
        append!(ys, y)
    end

    return (; ax, x=xs, y=ys)
end

function _gamma_shape_moment(df; occ::Float64=0.95)
    ratios = Float64[]

    for class in unique(df.class)
        sdf = df[df.class .== class, :]
        freqs = Meris.DataTools.get_frequencies(sdf; occ=occ, rescale=false)

        size(freqs, 2) == 0 && continue

        for col in eachcol(freqs)
            μ = mean(col)
            σ2 = var(col, corrected=false)

            (!isfinite(μ) || !isfinite(σ2) || μ <= 0) && continue

            push!(ratios, σ2 / (μ^2))
        end
    end

    isempty(ratios) && return nothing
    return mean(ratios)
end

function _add_icon!(parent_cell, icon_path;
    width=Relative(0.25), height=Relative(0.25),
    halign=0.08, valign=0.92
)
    axicon = Axis(
        parent_cell;
        width=width, height=height,
        halign=halign, valign=valign,
        tellwidth=false, tellheight=false,
    )

    icon = FileIO.load(icon_path)
    icon_small = imresize(icon, (256, 256))
    image!(axicon, rotr90(icon_small))
    hidedecorations!(axicon)
    hidespines!(axicon)

    return axicon
end

function _load_linguistic_df()
    df_arxiv = Meris.arXivLoader.load()
    df_arxiv.class .= "arx-" .* uppercase.(df_arxiv.class)
    select!(df_arxiv, :class, :sample_id, :component_id, :counts, :nreads)

    df_gut = Meris.GutenbergLoader.load()
    df_gut.class .= "guten-" .* uppercase.(df_gut.class)
    select!(df_gut, :class, :sample_id, :component_id, :counts, :nreads)

    df_rfc = Meris.RFCLoader.load()
    df_rfc.class .= uppercase.(df_rfc.class)
    select!(df_rfc, :class, :sample_id, :component_id, :counts, :nreads)

    return vcat(df_arxiv, df_gut, df_rfc)
end

function _load_microbial_df()
    df = Meris.OTULoader.load()
    select!(df, :class, :sample_id, :component_id, :counts, :nreads)
    return df
end

function _load_social_df()
    df_fin = Meris.FinanceLoader.load()
    df_fin = df_fin[endswith.(df_fin.class, "-daily"), :]
    df_fin.class = replace.(df_fin.class, "-daily" => "")
    df_fin.class .= "stock-" .* uppercase.(df_fin.class)
    select!(df_fin, :class, :sample_id, :component_id, :counts, :nreads)

    df_gow = Meris.GowallaLoader.load()
    df_gow.class .= "CHECK-IN"
    select!(df_gow, :class, :sample_id, :component_id, :counts, :nreads)

    df_lego = Meris.LEGOLoader.load()
    df_lego.class .= "LEGO"
    select!(df_lego, :class, :sample_id, :component_id, :counts, :nreads)

    return vcat(df_fin, df_gow, df_lego)
end

function _load_biology_df()
    df_bci = Meris.BCITreeLoader.load()
    df_bci.class .= "eco-BCI"
    select!(df_bci, :class, :sample_id, :component_id, :counts, :nreads)

    df_bio = Meris.BioTIMELoader.load()
    df_bio.class .= "eco-BT" .* string.(df_bio.class)
    select!(df_bio, :class, :sample_id, :component_id, :counts, :nreads)

    df_gtex = Meris.GTExLoader.load()
    df_gtex.class .= "gen-" .* string.(df_gtex.class)
    select!(df_gtex, :class, :sample_id, :component_id, :counts, :nreads)

    return vcat(df_gtex, df_bci, df_bio)
end

function _default_datasets()
    bases = [
        colorant"#1f77b4",
        colorant"#ff7f0e",
        colorant"#9467bd",
        colorant"#2ca02c",
        colorant"#d62728",
    ]

    return [
        (;
            key=:linguistic,
            title="Linguistic",
            loader=_load_linguistic_df,
            palette=shades(bases[1], 10),
            icon=joinpath(ICONDIR, "document.png"),
            icon_kw=(; width=Relative(0.18), height=Relative(0.18), halign=0.08, valign=0.92),
            occ=0.95,
            nbins=25,
        ),
        (;
            key=:microbial,
            title="Microbial",
            loader=_load_microbial_df,
            palette=shades(bases[2], 10),
            icon=joinpath(ICONDIR, "bacteria.png"),
            icon_kw=(; width=Relative(0.18), height=Relative(0.18), halign=0.08, valign=0.92),
            occ=0.95,
            nbins=25,
        ),
        (;
            key=:social,
            title="Social",
            loader=_load_social_df,
            palette=shades(bases[3], 8),
            icon=joinpath(ICONDIR, "socio-economic.png"),
            icon_kw=(; width=Relative(0.77 * 0.18), height=Relative(0.18), halign=0.08, valign=0.92),
            occ=0.95,
            nbins=25,
        ),
        (;
            key=:biology,
            title="Biology",
            loader=_load_biology_df,
            palette=vcat(shades(bases[4], 10)[1:7], shades(bases[5], 8)),
            icon=joinpath(ICONDIR, "eco.png"),
            icon_kw=(; width=Relative(0.18), height=Relative(0.18), halign=0.08, valign=0.92),
            occ=0.95,
            nbins=25,
        ),
    ]
end

function plot!(parent;
    datasets=_default_datasets(),
    font_scale::Float64=1.2,
    xlimits=(-12, 8),
    ylimits=(1e-5, 1.0),
    show_icons::Bool=true,
    panel_rowgap=5,
    panel_colgap=6,
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true])
    set_theme!(__theme)

    markers = [:circle, :rect, :diamond, :cross, :x, :utriangle, :dtriangle, :star4, :star6, :pentagon, :hexagon, :octagon]
    panel = GridLayout(parent)
    axes = Axis[]
    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]

    for (i, spec) in enumerate(datasets)
        ax = Axis(
            panel[positions[i]...];
            yscale=log10,
            limits=(xlimits[1], xlimits[2], ylimits[1], ylimits[2]),
            xlabelvisible=false,
            ylabelvisible=false,
            xticklabelsize=NATURE_TICK_PT * font_scale,
            yticklabelsize=NATURE_TICK_PT * font_scale,
            xgridvisible=false,
            ygridvisible=false,
        )

        df = spec.loader()
        nclasses = length(unique(df.class))
        colors = [spec.palette[mod1(j, length(spec.palette))] for j in 1:nclasses]
        afd_out = ax_afd(ax, df, colors, markers; nbins=spec.nbins, occ=spec.occ)

        xrange = -8:1e-2:5
        αmoment = _gamma_shape_moment(df; occ=spec.occ)
        if !isnothing(αmoment) && isfinite(αmoment)
            lines!(
                ax, xrange, exp.(Meris.LRDistr.lr_gamma(xrange, αmoment));
                color=:black,
                linestyle=:dash,
                linewidth=1.2,
            )
            text!(
                ax, 0.08, 0.72;
                space=:relative,
                text=L"\alpha = %$(round(αmoment, digits=2))",
                color=:black,
                fontsize=NATURE_TICK_PT * font_scale,
                align=(:left, :top),
            )
        end

        if show_icons && isfile(spec.icon)
            _add_icon!(panel[positions[i]...], spec.icon; spec.icon_kw...)
        end

        push!(axes, ax)
    end

    for ax in axes[2:end]
        linkxaxes!(axes[1], ax)
        linkyaxes!(axes[1], ax)
    end

    hidexdecorations!(axes[1]; grid=false)
    hidexdecorations!(axes[2]; grid=false)
    hideydecorations!(axes[2]; grid=false)
    hideydecorations!(axes[4]; grid=false)

    rowsize!(panel, 1, Relative(0.5))
    rowsize!(panel, 2, Relative(0.5))
    colsize!(panel, 1, Relative(0.5))
    colsize!(panel, 2, Relative(0.5))
    rowgap!(panel, panel_rowgap)
    colgap!(panel, panel_colgap)

    Label(
        panel[1:2, 0],
        L"p(z)";
        rotation=π / 2,
        fontsize=NATURE_AXIS_LABEL_PT * font_scale,
        tellheight=false,
    )
    Label(
        panel[3, 1:2],
        L"z";
        fontsize=NATURE_AXIS_LABEL_PT * font_scale,
        tellwidth=false,
    )

    return (; panel, axes)
end

function plot_afd(; ext="pdf", savefig::Bool=true, figname=nothing, kwargs...)
    fig = Figure(
        size=(0.62 * NATURE_DOUBLE_WIDTH_PT, 0.44 * NATURE_MAX_HEIGHT_PT),
        figure_padding=(6, 12, 10, 10),
    )

    plot!(fig[1, 1]; kwargs...)

    if savefig
        outfile = isnothing(figname) ? (Meris.FIGDIR * "afd.$ext") : figname
        save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

end # module AFDPlotter
