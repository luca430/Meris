#= Plot AFD-style panels from TL-prediction component-bin summaries. =#
module TLPredictionAFDPlotter

using Meris
using CairoMakie, MakiePublication, LaTeXStrings
using Colors
using DataFrames
using JLD2
using Random
using StatsBase

include("../scripts/module-scripts/macropatterns/TL_prediction.jl")
using .TLPrediction

const MM_TO_PT = 72.0 / 25.4
const NATURE_DOUBLE_WIDTH_PT = 183.0 * MM_TO_PT
const NATURE_MAX_HEIGHT_PT = 170.0 * MM_TO_PT
const NATURE_AXIS_LABEL_PT = 7
const NATURE_TICK_PT = 6
const NATURE_PANEL_LABEL_PT = 8
const NATURE_TEXT_PT = 6

function _shades(base::Colorant, n::Int)
    hsl = convert(HSL, base)
    return [HSL(hsl.h, hsl.s, l) for l in range(0.25, 0.85, length=n)]
end

function _default_datasets(; resultdir=Meris.DATADIR * "macro/tl-prediction/")
    bases = [
        colorant"#1f77b4",
        colorant"#ff7f0e",
        colorant"#9467bd",
        colorant"#2ca02c",
    ]

    return [
        (; key="linguistic", title="Linguistic", file=joinpath(resultdir, "linguistic.jld2"), palette=_shades(bases[1], 10)),
        (; key="microbial", title="Microbial", file=joinpath(resultdir, "microbial.jld2"), palette=_shades(bases[2], 10)),
        (; key="social", title="Social", file=joinpath(resultdir, "social.jld2"), palette=_shades(bases[3], 8)),
        (; key="biology", title="Biology", file=joinpath(resultdir, "biology.jld2"), palette=_shades(bases[4], 12)),
    ]
end

function _load_prediction_data(path::AbstractString)
    isfile(path) || error("Missing TL-prediction file: $path")
    d = JLD2.load(path)
    haskey(d, "component_bins") || error("File does not contain component_bins: $path")
    return d
end

function _downsampled_file(spec, prediction_data; downsampled_dir=Meris.DATADIR * "downsampled/")
    if haskey(prediction_data, "downsampled_file") && isfile(prediction_data["downsampled_file"])
        return prediction_data["downsampled_file"]
    end

    return joinpath(downsampled_dir, "$(spec.key).jld2")
end

function _load_downsampled_domain_df(spec, prediction_data)
    path = _downsampled_file(spec, prediction_data)
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

function _selected_bins_by_omega(component_bins::DataFrame, omega::Real; min_components::Int=10)
    omega > 0 || error("Omega must be positive, got $omega")
    omega_log = log10(float(omega))

    bin_counts = combine(
        groupby(component_bins, [:class, :coeff_bin, :bin_center_log, :bin_center]),
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

function _class_count_dfs(spec, prediction_data; downsample::Bool=false)
    occ = Float64(prediction_data["occ"])

    class_dfs = TLPrediction.divide_in_class(
        _load_downsampled_domain_df(spec, prediction_data);
        downsample=downsample,
        occ=occ,
    )

    out = Dict{String, DataFrame}()
    for (class, df) in class_dfs
        cdf = copy(df)
        cdf.class = fill(string(class), nrow(cdf))
        out[string(class)] = cdf
    end

    return out
end

function _afd_groups(component_bins::DataFrame, selected_rows::DataFrame, class_dfs::Dict{String, DataFrame};
    omega::Real,
    min_points::Int=10,
    minoccupancy::Float64=0.05,
)
    isempty(selected_rows.class) && return NamedTuple[]
    mean_threshold = omega^(-1)

    groups = NamedTuple[]
    markers = [:circle, :rect, :diamond, :cross, :x, :utriangle, :dtriangle, :star4, :star6, :pentagon, :hexagon, :octagon]

    for (i, row) in enumerate(eachrow(selected_rows))
        bin_components = component_bins[
            (component_bins.class .== row.class) .&
            (component_bins.coeff_bin .== row.coeff_bin),
            :,
        ]
        mean_col = :mean in propertynames(bin_components) ? :mean : :mean_B
        filter!(mean_col => m -> isfinite(m) && m > mean_threshold, bin_components)
        isempty(bin_components.component_id) && continue
        haskey(class_dfs, row.class) || continue

        keep = Set(bin_components.component_id)
        cdf = filter(r -> r.component_id in keep, class_dfs[row.class])
        nrow(cdf) == 0 && continue

        afd = Meris.AFD.compute(copy(cdf), :component_id; maxfrequency=1.0, minoccupancy=minoccupancy)
        length(afd.z) < min_points && continue

        push!(
            groups,
            (;
                z=afd.z,
                class=String(row.class),
                marker=markers[mod1(i, length(markers))],
            ),
        )
    end

    return groups
end

function _plot_afd_groups!(ax, groups, palette; nbins::Int=25)
    for (i, group) in enumerate(groups)
        x, y = Meris.DataTools.make_hist(group.z; nbins=nbins)
        scatter!(
            ax, x, y;
            color=:white,
            strokecolor=palette[mod1(i, length(palette))],
            marker=group.marker,
            markersize=6.5,
            strokewidth=0.75,
        )
    end

    return ax
end

function plot!(parent;
    datasets=_default_datasets(),
    omegas=(0.1, 1.0, 10.0),
    min_components::Int=10,
    min_points::Int=10,
    nbins::Int=25,
    downsample::Bool=false,
    font_scale::Float64=1.2,
    xlimits=(-4, 4),
    ylimits=(1e-3, 1.0),
)
    sc = Cycle([:color => :markercolor, :strokecolor => :color, :marker], covary=true)
    set_theme!(MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true, true]))

    panel = GridLayout(parent)
    nrows = length(omegas)
    ncols = length(datasets)

    prediction_data_by_key = Dict(spec.key => _load_prediction_data(spec.file) for spec in datasets)
    class_dfs_by_key = Dict(
        spec.key => _class_count_dfs(spec, prediction_data_by_key[spec.key]; downsample=downsample)
        for spec in datasets
    )

    for (r, omega) in enumerate(omegas)
        for (c, spec) in enumerate(datasets)
            ax = Axis(
                panel[r, c];
                yscale=log10,
                limits=(xlimits[1], xlimits[2], ylimits[1], ylimits[2]),
                xlabelvisible=r == nrows,
                ylabelvisible=c == 1,
                xticklabelsize=NATURE_TICK_PT * font_scale,
                yticklabelsize=NATURE_TICK_PT * font_scale,
                xgridvisible=false,
                ygridvisible=false,
            )

            prediction_data = prediction_data_by_key[spec.key]
            component_bins = prediction_data["component_bins"]
            minoccupancy = Float64(prediction_data["occ"])
            selected_rows = _selected_bins_by_omega(component_bins, omega; min_components=min_components)
            groups = _afd_groups(
                component_bins,
                selected_rows,
                class_dfs_by_key[spec.key];
                omega=omega,
                min_points=min_points,
                minoccupancy=minoccupancy,
            )
            _plot_afd_groups!(ax, groups, spec.palette; nbins=nbins)

            if r == 1
                Label(
                    panel[0, c],
                    spec.title;
                    fontsize=NATURE_AXIS_LABEL_PT * font_scale,
                    tellwidth=false,
                )
            end

            if c == 1
                Label(
                    panel[r, 0],
                    latexstring("\\Omega = ", string(omega));
                    rotation=π / 2,
                    fontsize=NATURE_TEXT_PT * font_scale,
                    tellheight=false,
                )
            end
        end
    end

    Label(
        panel[nrows + 1, 1:ncols],
        L"z";
        fontsize=NATURE_AXIS_LABEL_PT * font_scale,
        tellwidth=false,
    )
    Label(
        panel[1:nrows, -1],
        L"p(z)";
        rotation=π / 2,
        fontsize=NATURE_AXIS_LABEL_PT * font_scale,
        tellheight=false,
    )

    rowgap!(panel, 6)
    colgap!(panel, 8)

    return (; panel)
end

function plot_tl_prediction_afd(; ext="pdf", savefig::Bool=true, figname=nothing, kwargs...)
    fig = Figure(
        size=(1.2 * NATURE_DOUBLE_WIDTH_PT, 0.82 * NATURE_MAX_HEIGHT_PT),
        figure_padding=(8, 10, 8, 10),
    )

    plot!(fig[1, 1]; kwargs...)

    if savefig
        outfile = isnothing(figname) ? (Meris.FIGDIR * "tl-prediction-afd.$ext") : figname
        save(outfile, fig, pt_per_unit=1)
    end

    return fig
end

end # module TLPredictionAFDPlotter
