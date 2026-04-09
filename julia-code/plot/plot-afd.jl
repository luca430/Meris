#= Module to plot AFD panels in a Figure 2-style 2x2 layout =#
module AFDPlotter

using Meris
using DataFrames, StatsBase
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
const NATURE_PANEL_LABEL_PT = 8
const LINGUISTIC_STOPWORDS = Meris.arXivLoader.STOPWORDS

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

function _gamma_shape_moment(df; occ::Float64=0.999)
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

function _plot_gamma_fit!(ax, xrange, βmoment, label_fn;
    color=:black,
    linestyle=:dash,
    position=(0.08, 0.72),
    font_scale::Float64=1.0,
)
    if !isnothing(βmoment) && isfinite(βmoment)
        lines!(
            ax, xrange, exp.(Meris.LRDistr.lr_gamma(xrange, βmoment));
            color=color,
            linestyle=linestyle,
            linewidth=1.2,
        )
        text!(
            ax, position[1], position[2];
            space=:relative,
            text=label_fn(round(βmoment, digits=2)),
            color=color,
            fontsize=NATURE_TICK_PT * font_scale,
            align=(:left, :top),
        )
    end
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

    df = vcat(df_arxiv, df_gut, df_rfc)
    filter!(row -> !(lowercase(String(row.component_id)) in LINGUISTIC_STOPWORDS), df)
    return df
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

    df = vcat(df_gtex, df_bci, df_bio)
    sort!(df, :class; by=cls -> (startswith(cls, "gen-") ? 0 : 1, cls))
    return df
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
            occ=0.999,
            nbins=25,
        ),
        (;
            key=:microbial,
            title="Microbial",
            loader=_load_microbial_df,
            palette=shades(bases[2], 10),
            icon=joinpath(ICONDIR, "bacteria.png"),
            icon_kw=(; width=Relative(0.18), height=Relative(0.18), halign=0.08, valign=0.92),
            occ=0.999,
            nbins=25,
        ),
        (;
            key=:social,
            title="Social",
            loader=_load_social_df,
            palette=shades(bases[3], 8),
            icon=joinpath(ICONDIR, "socio-economic.png"),
            icon_kw=(; width=Relative(0.77 * 0.18), height=Relative(0.18), halign=0.08, valign=0.92),
            occ=0.999,
            nbins=25,
        ),
        (;
            key=:biology,
            title="Biology",
            loader=_load_biology_df,
            palette=vcat(shades(bases[5], 8)[1:4], shades(bases[4], 10)[1:7]),
            icon=joinpath(ICONDIR, "eco.png"),
            icon_kw=(; width=Relative(0.18), height=Relative(0.18), halign=0.08, valign=0.92),
            occ=0.999,
            nbins=25,
        ),
    ]
end

function _datasets_with_occ(datasets, occ::Float64)
    return [merge(spec, (; occ=occ)) for spec in datasets]
end

function _materialize_datasets(datasets)
    return [merge(spec, (; df=spec.loader())) for spec in datasets]
end

function _plot_subfigure!(parent;
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
    panel = parent isa GridLayout ? parent : GridLayout(parent)
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

        df = hasproperty(spec, :df) ? spec.df : spec.loader()
        nclasses = length(unique(df.class))
        colors = [spec.palette[mod1(j, length(spec.palette))] for j in 1:nclasses]
        afd_out = ax_afd(ax, df, colors, markers; nbins=spec.nbins, occ=spec.occ)

        xrange = -8:1e-2:5
        if spec.key == :biology
            df_gen = df[startswith.(df.class, "gen-"), :]
            df_eco = df[startswith.(df.class, "eco-"), :]
            _plot_gamma_fit!(ax, xrange, _gamma_shape_moment(df_gen; occ=spec.occ), β -> latexstring("\\beta_{\\mathrm{gen}} = ", string(β));
                color=:black,
                linestyle=:dash,
                position=(0.08, 0.72),
                font_scale=font_scale,
            )
            _plot_gamma_fit!(ax, xrange, _gamma_shape_moment(df_eco; occ=spec.occ), β -> latexstring("\\beta_{\\mathrm{eco}} = ", string(β));
                color=:black,
                linestyle=:dot,
                position=(0.08, 0.62),
                font_scale=font_scale,
            )
        else
            βmoment = _gamma_shape_moment(df; occ=spec.occ)
            _plot_gamma_fit!(ax, xrange, βmoment, β -> latexstring("\\beta = ", string(β));
                color=:black,
                linestyle=:dash,
                position=(0.08, 0.72),
                font_scale=font_scale,
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

function plot!(parent;
    datasets=_default_datasets(),
    occs=(0.999, 0.5),
    letters=('a', 'b'),
    font_scale::Float64=1.2,
    xlimits=(-12, 8),
    ylimits=(1e-5, 1.0),
    show_icons::Bool=true,
    panel_rowgap=5,
    panel_colgap=6,
    subfigure_colgap=12,
)
    container = GridLayout(parent)
    base_datasets = _materialize_datasets(datasets)

    for (i, occ) in enumerate(occs)
        subfig = GridLayout(container[1, i])
        _plot_subfigure!(subfig;
            datasets=_datasets_with_occ(base_datasets, occ),
            font_scale=font_scale,
            xlimits=xlimits,
            ylimits=ylimits,
            show_icons=show_icons,
            panel_rowgap=panel_rowgap,
            panel_colgap=panel_colgap,
        )
        if i <= length(letters)
            Label(
                container[1, i, TopLeft()],
                string(letters[i]);
                fontsize=NATURE_PANEL_LABEL_PT * font_scale,
                font=:bold,
                color=:black,
                halign=:left,
                valign=:bottom,
                padding=(0, 0, 6, 0),
            )
        end
    end

    colgap!(container, subfigure_colgap)
    return (; panel=container)
end

function plot_afd(; ext="pdf", savefig::Bool=true, figname=nothing, kwargs...)
    fig = Figure(
        size=(1.24 * NATURE_DOUBLE_WIDTH_PT, 0.44 * NATURE_MAX_HEIGHT_PT),
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
