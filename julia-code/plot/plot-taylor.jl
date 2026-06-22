#= Plot Taylor's law in log-log scale for all dataset classes. =#
module TaylorPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings
using JLD2
using DataFrames
using Statistics
using Random
using Colors, ColorTypes
using FileIO, ImageTransformations

import Meris

include("./colors/shadetester.jl")
using .Shades: shades

const ICONDIR = joinpath(Meris.FIGDIR, "icons")
const MM_TO_PT = 72.0 / 25.4
const NATURE_DOUBLE_WIDTH_PT = 183.0 * MM_TO_PT
const NATURE_MAX_HEIGHT_PT = 170.0 * MM_TO_PT
const NATURE_AXIS_LABEL_PT = 7
const NATURE_TICK_PT = 6
const NATURE_TEXT_PT = 6

const DOMAIN_BASES = (;
    linguistic=colorant"#1f77b4",
    microbial=colorant"#ff7f0e",
    social=colorant"#9467bd",
    ecological=colorant"#2ca02c",
    genetic=colorant"#d62728",
)

function _default_taylor_domains(TLDIR)
    return [
        (;
            key=:linguistic,
            title="Linguistic",
            file=joinpath(TLDIR, "linguistic.jld2"),
            class_filter=class -> true,
            palette=shades(DOMAIN_BASES.linguistic, 10),
            icon=joinpath(ICONDIR, "document.png"),
        ),
        (;
            key=:microbial,
            title="Microbial",
            file=joinpath(TLDIR, "microbial.jld2"),
            class_filter=class -> true,
            palette=shades(DOMAIN_BASES.microbial, 10),
            icon=joinpath(ICONDIR, "bacteria.png"),
        ),
        (;
            key=:social,
            title="Social",
            file=joinpath(TLDIR, "social.jld2"),
            class_filter=class -> true,
            palette=shades(DOMAIN_BASES.social, 8),
            icon=joinpath(ICONDIR, "socio-economic.png"),
        ),
        (;
            key=:ecological,
            title="Ecological",
            file=joinpath(TLDIR, "biology.jld2"),
            class_filter=class -> startswith(String(class), "eco-"),
            palette=shades(DOMAIN_BASES.ecological, 20),
            icon=joinpath(ICONDIR, "eco.png"),
        ),
        (;
            key=:genetic,
            title="Genetic",
            file=joinpath(TLDIR, "biology.jld2"),
            class_filter=class -> startswith(String(class), "gen-"),
            palette=shades(DOMAIN_BASES.genetic, 6),
            icon=joinpath(ICONDIR, "eco.png"),
        ),
    ]
end

function _default_downsampled_taylor_domains(; downsampled_dir=joinpath(Meris.DATADIR, "downsampled"))
    domains = _default_taylor_domains(joinpath(Meris.DATADIR, "macro", "taylor"))
    file_for_key = Dict(
        :linguistic => "linguistic.jld2",
        :microbial => "microbial.jld2",
        :social => "social.jld2",
        :ecological => "biology.jld2",
        :genetic => "biology.jld2",
    )

    return [
        merge(spec, (; downsampled_file=joinpath(downsampled_dir, file_for_key[spec.key])))
        for spec in domains
    ]
end

function _fig3_label_to_taylor_class(label)
    s = String(label)
    s = replace(s, "gutenberg-" => "guten-")
    s = replace(s, "GOWALLA" => "CHECK-IN")
    s = replace(s, "eco-BCI.Trees" => "eco-BCI")
    return s
end

function _fig3_classes_for_domain(domain::Symbol; FIG3DIR=joinpath(Meris.DATADIR, "fig3"))
    domain_dir = joinpath(FIG3DIR, domain == :ecological || domain == :genetic ? "biology" : String(domain))
    isdir(domain_dir) || return nothing

    classes = Set{String}()
    for file in sort(filter(f -> endswith(lowercase(f), ".jld2"), readdir(domain_dir)))
        out = JLD2.load(joinpath(domain_dir, file))["out"]
        for label in keys(out.pl)
            cls = _fig3_label_to_taylor_class(label)
            if domain == :ecological && !startswith(cls, "eco-")
                continue
            elseif domain == :genetic && !startswith(cls, "gen-")
                continue
            end
            push!(classes, cls)
        end
    end

    return classes
end

function _default_taylor_panels(domains)
    by_key = Dict(spec.key => spec for spec in domains)
    return [
        [by_key[:linguistic]],
        [by_key[:microbial]],
        [by_key[:social]],
        [by_key[:ecological], by_key[:genetic]],
    ]
end

function _load_tldf(path::AbstractString)
    isfile(path) || error("Missing Taylor data file: $path")
    data = JLD2.load(path)
    haskey(data, "tldf") || error("Taylor data file does not contain key `tldf`: $path")
    return data["tldf"]
end

function _compute_tldf_from_counts(df::DataFrame)
    dfs = DataFrame[]
    for class in sort(unique(df.class); by=string)
        sdf = df[df.class .== class, :]
        nrow(sdf) == 0 && continue
        count_df = copy(sdf)
        count_df.nreads .= one(eltype(count_df.nreads))
        tldf = Meris.Taylor.compute(count_df, :component_id; maxfrequency=Inf)
        tldf.class .= class
        push!(dfs, tldf)
    end

    isempty(dfs) && return DataFrame()
    return vcat(dfs...)
end

function _load_downsampled_tldf(spec; FIG3DIR=joinpath(Meris.DATADIR, "fig3"), use_fig3_classes::Bool=true)
    isfile(spec.downsampled_file) || error("Missing downsampled file: $(spec.downsampled_file)")
    data = JLD2.load(spec.downsampled_file)
    haskey(data, "ds_df") || error("Downsampled file does not contain `ds_df`: $(spec.downsampled_file)")

    df = data["ds_df"]
    select!(df, :class, :sample_id, :component_id, :counts, :nreads)
    df = filter(row -> spec.class_filter(row.class), df)

    if use_fig3_classes
        fig3_classes = _fig3_classes_for_domain(spec.key; FIG3DIR=FIG3DIR)
        if !isnothing(fig3_classes)
            df = filter(row -> String(row.class) in fig3_classes, df)
        end
    end

    return _compute_tldf_from_counts(df)
end

function _domain_df(spec; FIG3DIR=joinpath(Meris.DATADIR, "fig3"), use_fig3_classes::Bool=true)
    if hasproperty(spec, :downsampled_file)
        return _load_downsampled_tldf(spec; FIG3DIR=FIG3DIR, use_fig3_classes=use_fig3_classes)
    end

    df = _load_tldf(spec.file)
    df = filter(row -> spec.class_filter(row.class), df)

    if use_fig3_classes
        fig3_classes = _fig3_classes_for_domain(spec.key; FIG3DIR=FIG3DIR)
        if !isnothing(fig3_classes)
            df = filter(row -> String(row.class) in fig3_classes, df)
        end
    end

    return df
end

function _log_taylor_xy(df; mean_col=:omeanfrequency, var_col=:ovarfrequency)
    sdf = filter(row ->
        isfinite(row[mean_col]) &&
        isfinite(row[var_col]) &&
        row[mean_col] > 0 &&
        row[var_col] > 0,
        df
    )
    x, y = Meris.Taylor.clean_log(log10.(sdf[!, mean_col]), log10.(sdf[!, var_col]))
    return collect(x), collect(y)
end

function _log_taylor_series(df;
    bin_data::Bool=true,
    nbins::Int=25,
    min_points::Int=8,
    sample_fraction::Real=1.0,
    rng=Random.default_rng(),
    kwargs...
)
    x, y = _log_taylor_xy(df; kwargs...)
    length(x) < min_points && return Float64[], Float64[]
    if !bin_data && sample_fraction < 1.0
        n = max(min_points, ceil(Int, sample_fraction * length(x)))
        n = min(n, length(x))
        idx = randperm(rng, length(x))[1:n]
        x = x[idx]
        y = y[idx]
    end
    bin_data || return x, y

    nb = min(nbins, max(1, length(unique(x)) - 1))
    xb, yb = Meris.Taylor.binned_average(x, y; nbins=nb)
    xb, yb = Meris.Taylor.clean_log(xb, yb)
    return Float64.(xb), Float64.(yb)
end

function _linear_fit(x, y)
    length(x) >= 2 || return nothing
    mx = mean(x)
    my = mean(y)
    denom = sum(abs2, x .- mx)
    denom > 0 || return nothing
    slope = sum((x .- mx) .* (y .- my)) / denom
    intercept = my - slope * mx
    return (; intercept, slope)
end

function _omega_selected_series(df;
    omega::Real=1e-1,
    omega_half_width::Real=0.15,
    mean_col=:omeanfrequency,
    var_col=:ovarfrequency,
    sample_fraction::Real=1.0,
    min_points::Int=8,
    rng=Random.default_rng(),
)
    sdf = filter(row ->
        isfinite(row[mean_col]) &&
        isfinite(row[var_col]) &&
        row[mean_col] > 0 &&
        row[var_col] > row[mean_col],
        df
    )
    nrow(sdf) < min_points && return Float64[], Float64[]

    Ω = (sdf[!, var_col] .- sdf[!, mean_col]) ./ (sdf[!, mean_col] .^ 2)
    logΩ = log10.(Ω)
    target = log10(omega)
    mask = isfinite.(logΩ) .& (abs.(logΩ .- target) .<= omega_half_width)

    if count(mask) < min_points
        order = sortperm(abs.(logΩ .- target))
        n = min(max(min_points, count(isfinite.(logΩ))), length(order))
        mask .= false
        mask[order[1:n]] .= true
    end

    x = log10.(sdf[mask, mean_col])
    y = log10.(sdf[mask, var_col])
    x, y = Meris.Taylor.clean_log(x, y)

    if sample_fraction < 1.0 && length(x) > min_points
        n = max(min_points, ceil(Int, sample_fraction * length(x)))
        n = min(n, length(x))
        idx = randperm(rng, length(x))[1:n]
        x = x[idx]
        y = y[idx]
    end

    return Float64.(x), Float64.(y)
end

function _omega_curve(logxmin, logxmax, omega; n::Int=300)
    logxs = range(logxmin, logxmax; length=n)
    μ = exp10.(logxs)
    return logxs, log10.(omega .* μ .^ 2 .+ μ)
end

function _omega_label(omega::Real)
    e = log10(omega)
    if isapprox(e, round(e); atol=1e-8)
        return latexstring("\\Omega = 10^{", string(round(Int, e)), "}")
    end
    if isapprox(2e, round(2e); atol=1e-8)
        return latexstring("\\Omega = 10^{", string(round(e; digits=1)), "}")
    end
    return latexstring("\\Omega = ", string(round(omega; sigdigits=2)))
end

function _load_component_bins(path::AbstractString)
    isfile(path) || error("Missing TL-prediction file: $path")
    data = JLD2.load(path)
    haskey(data, "component_bins") || error("TL-prediction file does not contain `component_bins`: $path")
    df = copy(data["component_bins"])
    if :mean ∉ propertynames(df) && :mean_B in propertynames(df)
        df.mean = df.mean_B
    end
    if :var ∉ propertynames(df) && :var_B in propertynames(df)
        df.var = df.var_B
    end
    return df
end

function _default_omega_panels(; RESULTDIR=joinpath(Meris.DATADIR, "macro", "tl-prediction"))
    domains = _default_taylor_domains(joinpath(Meris.DATADIR, "macro", "taylor"))
    by_key = Dict(spec.key => spec for spec in domains)
    return [
        [(; by_key[:linguistic]..., prediction_file=joinpath(RESULTDIR, "linguistic.jld2"))],
        [(; by_key[:microbial]..., prediction_file=joinpath(RESULTDIR, "microbial.jld2"))],
        [(; by_key[:social]..., prediction_file=joinpath(RESULTDIR, "social.jld2"))],
        [
            (; by_key[:ecological]..., prediction_file=joinpath(RESULTDIR, "biology.jld2")),
            (; by_key[:genetic]..., prediction_file=joinpath(RESULTDIR, "biology.jld2")),
        ],
    ]
end

function _selected_omega_bin_rows(component_bins::DataFrame, classes, omega_log::Real; min_components::Int=10)
    d = filter(row ->
        String(row.class) in Set(String.(classes)) &&
        isfinite(row.mean) &&
        isfinite(row.var) &&
        row.mean > 0 &&
        row.var > 0 &&
        isfinite(row.bin_center_log),
        component_bins
    )
    isempty(d.class) && return DataFrame()

    bin_counts = combine(
        groupby(d, [:class, :coeff_bin, :bin_center_log]),
        nrow => :ncomponents,
    )
    filter!(:ncomponents => >=(min_components), bin_counts)
    isempty(bin_counts.class) && return bin_counts

    selected = combine(groupby(bin_counts, :class)) do df
        distances = abs.(df.bin_center_log .- omega_log)
        candidates = df[distances .== minimum(distances), :]
        candidates[argmax(candidates.ncomponents), :]
    end
    sort!(selected, :class)
    return selected
end

function _log_count_hist(counts; edges)
    c = Float64.(counts)
    c = c[isfinite.(c) .& (c .> 0)]
    length(c) == 0 && return Float64[]

    bin_idx = searchsortedlast.(Ref(edges), c)
    bin_idx = clamp.(bin_idx, 1, length(edges) - 1)
    y = zeros(length(edges) - 1)
    for idx in bin_idx
        y[idx] += 1
    end
    widths = diff(edges)
    y ./= (sum(y) .* widths)
    return y
end

function _gutenberg_en_color()
    classes = ["RFC", "arx-ASTRO-PH", "arx-MATH", "arx-PHYSICS", "arx-Q-BIO", "guten-EN", "guten-IT"]
    palette = shades(DOMAIN_BASES.linguistic, 10)
    idx = findfirst(==("guten-EN"), classes)
    return palette[mod1(idx, length(palette))]
end

function _class_color(spec, class_name, classes)
    idx = findfirst(==(String(class_name)), String.(classes))
    isnothing(idx) && error("Class $class_name not found for $(spec.key)")
    return spec.palette[mod1(idx, length(spec.palette))]
end

function _single_count_distribution_specs(; downsampled_dir=joinpath(Meris.DATADIR, "downsampled"))
    domains = _default_downsampled_taylor_domains(; downsampled_dir=downsampled_dir)
    by_key = Dict(spec.key => spec for spec in domains)
    return [
        merge(by_key[:linguistic], (; class_name="guten-EN", panel_title="guten-EN")),
        merge(by_key[:microbial], (; class_name="GUT1", panel_title="GUT1")),
        merge(by_key[:social], (; class_name="CHECK-IN", panel_title="CHECK-IN")),
        merge(by_key[:ecological], (; class_name="eco-BCI", panel_title="eco-BCI")),
        merge(by_key[:genetic], (; class_name="gen-BRAIN", panel_title="gen-BRAIN")),
    ]
end

function _taylor_markers()
    return [:circle, :rect, :diamond, :utriangle, :dtriangle, :cross, :x, :star4, :pentagon, :hexagon]
end

function _class_order(classes)
    return sort(collect(classes); by=class -> (
        startswith(String(class), "gen-") ? 0 :
        startswith(String(class), "eco-") ? 1 : 2,
        String(class),
    ))
end

function _add_icon!(parent_cell, icon_path;
    width=Relative(0.17), height=Relative(0.17),
    halign=0.97, valign=0.03
)
    isfile(icon_path) || return nothing

    axicon = Axis(
        parent_cell;
        width=width,
        height=height,
        halign=halign,
        valign=valign,
        aspect=DataAspect(),
        tellwidth=false,
        tellheight=false,
    )

    icon = FileIO.load(icon_path)
    icon_small = imresize(icon, (256, 256))
    image!(axicon, rotr90(icon_small))
    hidedecorations!(axicon)
    hidespines!(axicon)

    return axicon
end

function _fit_label(fit)
    return latexstring("b = ", string(round(fit.slope; digits=2)))
end

"""
    plot!(parent; kwargs...)

Plot Taylor's law for every dataset class in the Taylor macro data. Points are
bin averages in log10(mean abundance) on the x-axis and log10(variance) on the
y-axis. Each dataset class receives its own linear fit, while colors encode the
five requested domains: linguistic, microbial, social, ecological, and genetic.
"""
function plot!(parent;
    TLDIR=joinpath(Meris.DATADIR, "macro", "taylor"),
    FIG3DIR=joinpath(Meris.DATADIR, "fig3"),
    domains=_default_taylor_domains(TLDIR),
    panels=_default_taylor_panels(domains),
    use_fig3_classes::Bool=true,
    bin_data::Bool=true,
    nbins::Int=25,
    sample_fraction::Real=1.0,
    seed::Int=1234,
    min_points::Int=8,
    mean_col=:omeanfrequency,
    var_col=:ovarfrequency,
    font_scale::Float64=1.7,
    markersize::Real=7.5,
    fit_linewidth::Real=1.15,
    vertical_shift_step::Real=11.0,
    show_icons::Bool=true,
    show_fits::Bool=true,
    show_fit_labels::Bool=true,
    show_reference_lines::Bool=false,
    reference_linewidth::Real=2.0,
    reference_anchor_shift=(-0.9, -0.75),
    reference_gray_shift::Real=-0.35,
    reference_black_shift::Real=0.35,
    axis_ylimits=nothing,
    limits=(nothing, nothing, nothing, nothing),
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    set_theme!(MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true]))

    markers = _taylor_markers()
    panel = GridLayout(parent)
    axes = Axis[]
    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]
    plotted_x = Float64[]
    rng = MersenneTwister(seed)
    df_cache = Dict{Symbol, DataFrame}()
    get_domain_df(spec) = get!(df_cache, spec.key) do
        _domain_df(spec; FIG3DIR=FIG3DIR, use_fig3_classes=use_fig3_classes)
    end

    for (i, panel_specs) in enumerate(panels)
        ax = Axis(
            panel[positions[i]...];
            xlabel=L"\log_{10}\,\mu",
            ylabel=L"\log_{10}\,\sigma^2",
            xlabelsize=NATURE_AXIS_LABEL_PT * font_scale,
            ylabelsize=NATURE_AXIS_LABEL_PT * font_scale,
            xticklabelsize=NATURE_TICK_PT * font_scale,
            yticklabelsize=NATURE_TICK_PT * font_scale,
            limits=limits,
        )

        fits = NamedTuple[]
        panel_x = Float64[]
        panel_y = Float64[]
        nseries = sum(length(unique(get_domain_df(spec).class)) for spec in panel_specs)
        series_index = 0
        for spec in panel_specs
            df = get_domain_df(spec)
            classes = _class_order(unique(df.class))
            colors = [spec.palette[mod1(j, length(spec.palette))] for j in eachindex(classes)]

            for (j, class) in enumerate(classes)
                sdf = df[df.class .== class, :]
                x, y = _log_taylor_series(
                    sdf;
                    bin_data=bin_data,
                    nbins=nbins,
                    sample_fraction=sample_fraction,
                    rng=rng,
                    min_points=min_points,
                    mean_col=mean_col,
                    var_col=var_col,
                )
                length(x) < 2 && continue
                append!(plotted_x, x)

                series_index += 1
                # Additive in log10 space; multiplicative in the original variance scale.
                yshift = (series_index - (nseries + 1) / 2) * vertical_shift_step
                ydisplay = y .+ yshift
                append!(panel_x, x)
                append!(panel_y, ydisplay)
                color = colors[j]
                marker = markers[mod1(j, length(markers))]
                scatter!(
                    ax,
                    x,
                    ydisplay;
                    color=:white,
                    strokecolor=color,
                    marker=marker,
                    markersize=markersize,
                    strokewidth=0.55,
                )

                fit = _linear_fit(x, y)
                isnothing(fit) && continue
                display_fit = (; intercept=fit.intercept + yshift, slope=fit.slope)

                xs = range(minimum(x), maximum(x); length=80)
                if show_fits
                    lines!(
                        ax,
                        xs,
                        display_fit.intercept .+ display_fit.slope .* xs;
                        color=:black,
                        linewidth=fit_linewidth,
                        linestyle=(:dash, :dense),
                    )
                end
                right_idx = argmax(x)
                xpad = 0.08 * (maximum(x) - minimum(x))
                label_x = x[right_idx] + xpad
                push!(plotted_x, label_x)
                push!(fits, (; class, fit, label_x, label_y=ydisplay[right_idx]))
            end

            GC.gc()
        end

        if show_reference_lines && !isempty(panel_x) && !isempty(panel_y)
            x0 = minimum(panel_x) + reference_anchor_shift[1]
            y0 = minimum(panel_y) + reference_anchor_shift[2]
            xcross = x0 + reference_gray_shift - reference_black_shift
            xleft = min(x0, xcross) - 0.05 * (maximum(panel_x) - minimum(panel_x))
            xs = range(xleft, maximum(panel_x); length=140)
            lines!(
                ax,
                xs,
                y0 .+ reference_gray_shift .+ (xs .- x0);
                color=:gray45,
                linewidth=reference_linewidth,
                linestyle=(:dash, :dense),
            )
            text!(
                ax,
                last(xs),
                y0 + reference_gray_shift + (last(xs) - x0);
                text=L"b = 1",
                align=(:left, :center),
                fontsize=NATURE_TEXT_PT * font_scale,
                color=:gray45,
            )
            lines!(
                ax,
                xs,
                y0 .+ reference_black_shift .+ 2 .* (xs .- x0);
                color=:black,
                linewidth=reference_linewidth,
                linestyle=(:dash, :dense),
            )
            text!(
                ax,
                last(xs),
                y0 + reference_black_shift + 2 * (last(xs) - x0);
                text=L"b = 2",
                align=(:left, :center),
                fontsize=NATURE_TEXT_PT * font_scale,
                color=:black,
            )
        end

        if show_fit_labels
            for f in fits
                text!(
                    ax,
                    f.label_x,
                    f.label_y;
                    text=_fit_label(f.fit),
                    align=(:left, :center),
                    fontsize=NATURE_TEXT_PT * font_scale * 1.05,
                    color=:black,
                )
            end
        end

        if show_icons
            icon_scale = any(spec.key in (:ecological, :genetic) for spec in panel_specs) ? 0.21 : 0.17
            _add_icon!(
                panel[positions[i]...],
                panel_specs[1].icon;
                width=Relative(icon_scale),
                height=Relative(icon_scale),
            )
        end

        push!(axes, ax)
    end

    for ax in axes[2:end]
        linkxaxes!(axes[1], ax)
        if isnothing(axis_ylimits)
            linkyaxes!(axes[1], ax)
        end
    end

    if !isnothing(axis_ylimits)
        for (i, ylim) in enumerate(axis_ylimits)
            i > length(axes) && break
            isnothing(ylim) && continue
            ylims!(axes[i], ylim...)
        end
    end

    if !isempty(plotted_x) && isnothing(limits[1]) && isnothing(limits[2])
        xmin, xmax = extrema(plotted_x)
        xspan = xmax - xmin
        xlims!(axes[1], xmin - 0.05 * xspan, xmax + 0.34 * xspan)
    end

    hideydecorations!(axes[2]; grid=false)
    hideydecorations!(axes[4]; grid=false)
    hidexdecorations!(axes[1]; grid=false)
    hidexdecorations!(axes[2]; grid=false)

    colgap!(panel, 8)
    rowgap!(panel, 8)

    return (; panel, axes)
end

function plot_taylor(;
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    kwargs...
)
    fig = Figure(
        size=(NATURE_DOUBLE_WIDTH_PT, 0.78 * NATURE_MAX_HEIGHT_PT),
        figure_padding=(10, 10, 12, 12),
    )

    plot!(fig[1, 1]; kwargs...)

    if savefig
        outfile = isnothing(figname) ? joinpath(Meris.FIGDIR, "taylor.$ext") : figname
        save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

function plot_taylor_unbinned(;
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    kwargs...
)
    outfile = isnothing(figname) ? joinpath(Meris.FIGDIR, "taylor-unbinned.$ext") : figname

    return plot_taylor(;
        ext=ext,
        savefig=savefig,
        figname=outfile,
        bin_data=false,
        sample_fraction=1 / 3,
        kwargs...
    )
end

function plot_taylor_reference(;
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    kwargs...
)
    outfile = isnothing(figname) ? joinpath(Meris.FIGDIR, "taylor-reference.$ext") : figname

    return plot_taylor(;
        ext=ext,
        savefig=savefig,
        figname=outfile,
        vertical_shift_step=0.0,
        show_fits=false,
        show_fit_labels=false,
        show_reference_lines=true,
        axis_ylimits=[(-16, 2), (-16, 2), nothing, nothing],
        kwargs...
    )
end

function _downsampled_plot_args(; downsampled_dir=joinpath(Meris.DATADIR, "downsampled"))
    domains = _default_downsampled_taylor_domains(; downsampled_dir=downsampled_dir)
    return (; domains, panels=_default_taylor_panels(domains))
end

function plot_taylor_downsampled(;
    downsampled_dir=joinpath(Meris.DATADIR, "downsampled"),
    kwargs...
)
    return plot_taylor(; _downsampled_plot_args(; downsampled_dir=downsampled_dir)..., kwargs...)
end

function plot_taylor_unbinned_downsampled(;
    downsampled_dir=joinpath(Meris.DATADIR, "downsampled"),
    kwargs...
)
    return plot_taylor_unbinned(; _downsampled_plot_args(; downsampled_dir=downsampled_dir)..., kwargs...)
end

function plot_taylor_reference_downsampled(;
    downsampled_dir=joinpath(Meris.DATADIR, "downsampled"),
    axis_ylimits=nothing,
    kwargs...
)
    return plot_taylor_reference(;
        _downsampled_plot_args(; downsampled_dir=downsampled_dir)...,
        axis_ylimits=axis_ylimits,
        kwargs...
    )
end

function plot_taylor_omega(;
    omega=nothing,
    omegas=10.0 .^ [-1.0, -0.5, 0.0, 0.5, 1.0],
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    RESULTDIR=joinpath(Meris.DATADIR, "macro", "tl-prediction"),
    FIG3DIR=joinpath(Meris.DATADIR, "fig3"),
    panels=_default_omega_panels(; RESULTDIR=RESULTDIR),
    use_fig3_classes::Bool=true,
    bin_data::Bool=false,
    nbins::Int=25,
    sample_fraction::Real=1 / 3,
    seed::Int=1234,
    min_points::Int=8,
    font_scale::Float64=1.7,
    markersize::Real=7.5,
    strokewidth::Real=0.55,
    curve_linewidth::Real=1.15,
    omega_shift_step::Real=5.5,
    omega_label_xpad_fraction::Real=0.03,
    omega_label_font_scale::Real=0.90,
    show_icons::Bool=true,
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    set_theme!(MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true]))

    markers = _taylor_markers()
    omega_values = isnothing(omega) ? collect(omegas) : [omega]
    omega_shifts = [
        (k - (length(omega_values) + 1) / 2) * omega_shift_step
        for k in eachindex(omega_values)
    ]
    fig = Figure(
        size=(NATURE_DOUBLE_WIDTH_PT, 0.78 * NATURE_MAX_HEIGHT_PT),
        figure_padding=(10, 10, 12, 12),
    )
    grid = GridLayout(fig[1, 1])
    axes = Axis[]
    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]
    rng = MersenneTwister(seed)
    plotted_x = Float64[]

    for (i, panel_specs) in enumerate(panels)
        ax = Axis(
            grid[positions[i]...];
            xlabel=L"\log_{10}\,\mu",
            ylabel=L"\log_{10}\,\sigma^2",
            xlabelsize=NATURE_AXIS_LABEL_PT * font_scale,
            ylabelsize=NATURE_AXIS_LABEL_PT * font_scale,
            xticklabelsize=NATURE_TICK_PT * font_scale,
            yticklabelsize=NATURE_TICK_PT * font_scale,
        )

        panel_x = Float64[]
        panel_y = Float64[]
        for spec in panel_specs
            component_bins = _load_component_bins(spec.prediction_file)
            classes = _class_order(unique(component_bins.class))
            if use_fig3_classes
                fig3_classes = _fig3_classes_for_domain(spec.key; FIG3DIR=FIG3DIR)
                if !isnothing(fig3_classes)
                    classes = [class for class in classes if String(class) in fig3_classes]
                end
            end
            classes = [class for class in classes if spec.class_filter(class)]
            class_colors = [spec.palette[mod1(j, length(spec.palette))] for j in eachindex(classes)]
            for (k, Ω) in enumerate(omega_values)
                selected_rows = _selected_omega_bin_rows(component_bins, classes, log10(Ω); min_components=min_points)

                for row in eachrow(selected_rows)
                    d = component_bins[
                        (component_bins.class .== row.class) .&
                        (component_bins.coeff_bin .== row.coeff_bin),
                        :,
                    ]
                    filter!([:mean, :var] => (m, v) -> isfinite(m) && isfinite(v) && m > 0 && v > 0, d)
                    nrow(d) < min_points && continue

                    if sample_fraction < 1.0
                        n = max(min_points, ceil(Int, sample_fraction * nrow(d)))
                        n = min(n, nrow(d))
                        d = d[randperm(rng, nrow(d))[1:n], :]
                    end

                    x = log10.(d.mean)
                    y = log10.(d.var) .+ omega_shifts[k]
                    x, y = Meris.Taylor.clean_log(x, y)
                    if bin_data
                        nb = min(nbins, max(1, length(unique(x)) - 1))
                        x, y = Meris.Taylor.binned_average(x, y; nbins=nb)
                        x, y = Meris.Taylor.clean_log(x, y)
                    end
                    class_index = findfirst(==(row.class), classes)
                    isnothing(class_index) && continue

                    append!(panel_x, x)
                    append!(panel_y, y)
                    append!(plotted_x, x)
                    scatter!(
                        ax,
                        x,
                        y;
                        color=:white,
                        strokecolor=class_colors[class_index],
                        marker=markers[mod1(class_index, length(markers))],
                        markersize=markersize,
                        strokewidth=strokewidth,
                    )
                end
            end

            component_bins = nothing
            GC.gc()
        end

        if !isempty(panel_x)
            for (k, Ω) in enumerate(omega_values)
                logxs, logys = _omega_curve(minimum(panel_x), maximum(panel_x), Ω)
                logys = logys .+ omega_shifts[k]
                lines!(
                    ax,
                    logxs,
                    logys;
                    color=:black,
                    linewidth=curve_linewidth,
                    linestyle=(:dash, :dense),
                )
                text!(
                    ax,
                    last(logxs) + omega_label_xpad_fraction * (maximum(panel_x) - minimum(panel_x)),
                    last(logys);
                    text=_omega_label(Ω),
                    align=(:left, :center),
                    fontsize=NATURE_TEXT_PT * font_scale * omega_label_font_scale,
                    color=:black,
                )
            end
        end

        if show_icons
            icon_scale = any(spec.key in (:ecological, :genetic) for spec in panel_specs) ? 0.21 : 0.17
            _add_icon!(
                grid[positions[i]...],
                panel_specs[1].icon;
                width=Relative(icon_scale),
                height=Relative(icon_scale),
            )
        end

        push!(axes, ax)
    end

    for ax in axes[2:end]
        linkxaxes!(axes[1], ax)
    end
    if !isempty(plotted_x)
        xmin, xmax = extrema(plotted_x)
        xspan = xmax - xmin
        xlims!(axes[1], xmin - 0.05 * xspan, xmax + 0.32 * xspan)
    end

    hideydecorations!(axes[2]; grid=false)
    hideydecorations!(axes[4]; grid=false)
    hidexdecorations!(axes[1]; grid=false)
    hidexdecorations!(axes[2]; grid=false)
    colgap!(grid, 8)
    rowgap!(grid, 8)

    if savefig
        outfile = isnothing(figname) ? joinpath(Meris.FIGDIR, "taylor-omega.$ext") : figname
        save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

function plot_taylor_omega_binned(;
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    kwargs...
)
    outfile = isnothing(figname) ? joinpath(Meris.FIGDIR, "taylor-omega-binned.$ext") : figname

    return plot_taylor_omega(;
        ext=ext,
        savefig=savefig,
        figname=outfile,
        bin_data=true,
        sample_fraction=1.0,
        kwargs...
    )
end

function plot_gutenberg_en_count_distribution(;
    downsampled_file=joinpath(Meris.DATADIR, "downsampled", "linguistic.jld2"),
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    class_name="guten-EN",
    nbins::Int=55,
    font_scale::Float64=1.8,
)
    data = JLD2.load(downsampled_file)
    haskey(data, "ds_df") || error("Downsampled file does not contain `ds_df`: $downsampled_file")
    df = data["ds_df"]
    df = df[df.class .== class_name, :]
    nrow(df) > 0 || error("No rows found for class $class_name in $downsampled_file")

    counts = Float64.(df.counts)
    cmin = minimum(counts[counts .> 0])
    cmax = maximum(counts)
    edges = exp10.(range(log10(cmin), log10(cmax); length=nbins + 1))
    centers = sqrt.(edges[1:end-1] .* edges[2:end])

    fig = Figure(
        size=(0.82 * NATURE_DOUBLE_WIDTH_PT, 0.62 * NATURE_MAX_HEIGHT_PT),
        figure_padding=(12, 12, 12, 12),
    )
    ax = Axis(
        fig[1, 1];
        xlabel=L"\text{count } n",
        ylabel=L"p(n)",
        xscale=log10,
        yscale=log10,
        xlabelsize=NATURE_AXIS_LABEL_PT * font_scale,
        ylabelsize=NATURE_AXIS_LABEL_PT * font_scale,
        xticklabelsize=NATURE_TICK_PT * font_scale,
        yticklabelsize=NATURE_TICK_PT * font_scale,
    )

    hists = Vector{Float64}[]
    for sdf in groupby(df, :sample_id)
        y = _log_count_hist(sdf.counts; edges=edges)
        isempty(y) && continue
        push!(hists, y)
        mask = y .> 0
        lines!(
            ax,
            centers[mask],
            y[mask];
            color=(:gray50, 0.38),
            linewidth=0.28,
        )
    end

    isempty(hists) && error("No sample histograms could be computed for $class_name")
    avg_y = mean(reduce(hcat, hists); dims=2)[:, 1]
    mask = avg_y .> 0
    lines!(
        ax,
        centers[mask],
        avg_y[mask];
        color=_gutenberg_en_color(),
        linewidth=2.4,
    )

    if savefig
        outfile = isnothing(figname) ? joinpath(Meris.FIGDIR, "gutenberg-en-count-distribution.$ext") : figname
        save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

function _plot_count_distribution_panel!(ax, df, class_name, color;
    nbins::Int=55,
    thin_color=(:gray50, 0.38),
    thin_linewidth::Real=0.28,
    avg_linewidth::Real=2.4,
)
    sdf = df[df.class .== class_name, :]
    nrow(sdf) > 0 || error("No rows found for class $class_name")

    counts = Float64.(sdf.counts)
    positive_counts = counts[counts .> 0]
    isempty(positive_counts) && error("No positive counts found for class $class_name")

    cmin = minimum(positive_counts)
    cmax = maximum(positive_counts)
    edges = exp10.(range(log10(cmin), log10(cmax); length=nbins + 1))
    centers = sqrt.(edges[1:end-1] .* edges[2:end])

    hists = Vector{Float64}[]
    for sample_df in groupby(sdf, :sample_id)
        y = _log_count_hist(sample_df.counts; edges=edges)
        isempty(y) && continue
        push!(hists, y)
        mask = y .> 0
        lines!(
            ax,
            centers[mask],
            y[mask];
            color=thin_color,
            linewidth=thin_linewidth,
        )
    end

    isempty(hists) && error("No sample histograms could be computed for $class_name")
    avg_y = mean(reduce(hcat, hists); dims=2)[:, 1]
    mask = avg_y .> 0
    lines!(
        ax,
        centers[mask],
        avg_y[mask];
        color=color,
        linewidth=avg_linewidth,
    )

    return nothing
end

function plot_single_domain_count_distributions(;
    downsampled_dir=joinpath(Meris.DATADIR, "downsampled"),
    specs=_single_count_distribution_specs(; downsampled_dir=downsampled_dir),
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    nbins::Int=55,
    font_scale::Float64=1.55,
)
    set_theme!(MakiePublication.theme_acs())

    fig = Figure(
        size=(NATURE_DOUBLE_WIDTH_PT, 0.95 * NATURE_MAX_HEIGHT_PT),
        figure_padding=(10, 10, 12, 12),
    )
    grid = GridLayout(fig[1, 1])
    positions = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2)]
    axes = Axis[]

    for (i, spec) in enumerate(specs)
        ax = Axis(
            grid[positions[i]...];
            title=spec.panel_title,
            xlabel=L"\text{count } n",
            ylabel=L"p(n)",
            xscale=log10,
            yscale=log10,
            titlesize=NATURE_TEXT_PT * font_scale,
            xlabelsize=NATURE_AXIS_LABEL_PT * font_scale,
            ylabelsize=NATURE_AXIS_LABEL_PT * font_scale,
            xticklabelsize=NATURE_TICK_PT * font_scale,
            yticklabelsize=NATURE_TICK_PT * font_scale,
        )

        data = JLD2.load(spec.downsampled_file)
        haskey(data, "ds_df") || error("Downsampled file does not contain `ds_df`: $(spec.downsampled_file)")
        df = data["ds_df"]
        df = filter(row -> spec.class_filter(row.class), df)
        classes = _class_order(unique(df.class))
        color = _class_color(spec, spec.class_name, classes)
        _plot_count_distribution_panel!(ax, df, spec.class_name, color; nbins=nbins)
        push!(axes, ax)
    end

    hideydecorations!(axes[2]; grid=false)
    hideydecorations!(axes[3]; grid=false)
    hideydecorations!(axes[5]; grid=false)

    colgap!(grid, 8)
    rowgap!(grid, 8)

    if savefig
        outfile = isnothing(figname) ? joinpath(Meris.FIGDIR, "single-domain-count-distributions.$ext") : figname
        save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

function _dataset_display_name(class_name)
    s = String(class_name)
    s = replace(s, "arx-" => "arXiv ")
    s = replace(s, "guten-" => "Gutenberg ")
    s = replace(s, "stock-" => "Stock ")
    s = replace(s, "eco-" => "")
    s = replace(s, "gen-" => "GTEx ")
    return s
end

function _dataset_entries(spec; FIG3DIR=joinpath(Meris.DATADIR, "fig3"))
    classes = _fig3_classes_for_domain(spec.key; FIG3DIR=FIG3DIR)
    isnothing(classes) && (classes = String[])
    classes = [class for class in _class_order(classes) if spec.class_filter(class)]
    markers = _taylor_markers()

    return [
        (;
            class=class,
            label=_dataset_display_name(class),
            color=spec.palette[mod1(j, length(spec.palette))],
            marker=markers[mod1(j, length(markers))],
        )
        for (j, class) in enumerate(classes)
    ]
end

function _plot_dataset_key_panel!(parent_cell, panel_specs;
    FIG3DIR=joinpath(Meris.DATADIR, "fig3"),
    font_scale::Float64=1.7,
    markersize::Real=13,
)
    entries = reduce(vcat, [_dataset_entries(spec; FIG3DIR=FIG3DIR) for spec in panel_specs])
    n = length(entries)
    ax = Axis(
        parent_cell;
        limits=(0, 1, 0, max(n + 1, 2)),
        xautolimitmargin=(0, 0),
        yautolimitmargin=(0, 0),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    for (i, entry) in enumerate(entries)
        y = n - i + 1
        scatter!(
            ax,
            [0.13],
            [y];
            color=:white,
            strokecolor=entry.color,
            marker=entry.marker,
            markersize=markersize,
            strokewidth=0.85,
        )
        text!(
            ax,
            0.22,
            y;
            text=entry.label,
            align=(:left, :center),
            fontsize=NATURE_TEXT_PT * font_scale,
            color=:black,
        )
    end

    return ax
end

function plot_dataset_key(;
    TLDIR=joinpath(Meris.DATADIR, "macro", "taylor"),
    FIG3DIR=joinpath(Meris.DATADIR, "fig3"),
    domains=_default_taylor_domains(TLDIR),
    panels=_default_taylor_panels(domains),
    ext="pdf",
    savefig::Bool=true,
    figname=nothing,
    font_scale::Float64=1.75,
    markersize::Real=13,
)
    set_theme!(MakiePublication.theme_acs())

    fig = Figure(
        size=(NATURE_DOUBLE_WIDTH_PT, 0.62 * NATURE_MAX_HEIGHT_PT),
        figure_padding=(10, 10, 10, 10),
    )
    grid = GridLayout(fig[1, 1])
    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]

    for (i, panel_specs) in enumerate(panels)
        _plot_dataset_key_panel!(
            grid[positions[i]...],
            panel_specs;
            FIG3DIR=FIG3DIR,
            font_scale=font_scale,
            markersize=markersize,
        )
    end

    colgap!(grid, 8)
    rowgap!(grid, 8)

    if savefig
        outfile = isnothing(figname) ? joinpath(Meris.FIGDIR, "dataset-key.$ext") : figname
        save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

end # module TaylorPlotter

if abspath(PROGRAM_FILE) == @__FILE__
    using Meris
    using .TaylorPlotter

    TaylorPlotter.plot_taylor()
end
