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
function prepare(;
        categories=["linguistic", "microbial", "social", "biology"]
    )

    ## LINGUISTIC ##
    if "linguistic" in categories
        
        @info "Loading linguistic data..."
        
        #| arXiv |#
        df_arxiv = Meris.arXivLoader.load()
        df_arxiv.class .= "arx-" .* uppercase.(df_arxiv.class)
        select!(df_arxiv, :class, :sample_id, :component_id, :counts, :nreads)
        
        #| Gutenberg |#
        df_gut = Meris.GutenbergLoader.load()
        df_gut.class = "guten-" .* uppercase.(df_gut.class)
        select!(df_gut, :class, :sample_id, :component_id, :counts, :nreads)
        
        #| RFCs |#
        df_rfc = Meris.RFCLoader.load()
        df_rfc.class .= uppercase.(df_rfc.class)
        select!(df_rfc, :class, :sample_id, :component_id, :counts, :nreads)
        
        df = vcat(df_arxiv, df_gut, df_rfc)
        @info "Working linguistic data..."
        Taylor.compute(df; save=true, filename=Meris.DATADIR * "macro/taylor/linguistic.jld2")
    
        df = nothing
        df_arxiv = df_gut = df_rfc = nothing
        GC.gc()
    
    ## MICROBIAL ##
    elseif "microbial" in categories
        
        @info "Loading microbial data..."
        
        #| OTU |#
        df_otu = Meris.OTULoader.load()
        select!(df_otu, :class, :sample_id, :component_id, :counts, :nreads)
        
        df = vcat(df_otu)
        @info "Working microbial data..."
       Taylor.compute(df; save=true, filename=Meris.DATADIR * "macro/taylor/microbial.jld2")
    
        df = nothing
        df_otu = nothing
        GC.gc()
    
    ## SOCIAL ##
    elseif "social" in categories
        
        @info "Loading social data..."
        
        #| FINANCE |#
        df_fin = Meris.FinanceLoader.load()
        df_fin = df_fin[endswith.(df_fin.class, "-daily"), :]
        df_fin.class = replace.(df_fin.class, "-daily" => "")
        df_fin.class .= "stock-" .* uppercase.(df_fin.class)
        select!(df_fin, :class, :sample_id, :component_id, :counts, :nreads)
        
        #| Gowalla |#
        df_gow = Meris.GowallaLoader.load()
        df_gow.class .= "CHECK-IN"
        select!(df_gow, :class, :sample_id, :component_id, :counts, :nreads)
        
        #| LEGO |#
        df_lego = Meris.LEGOLoader.load()
        df_lego.class .= "LEGO"
        select!(df_lego, :class, :sample_id, :component_id, :counts, :nreads)
        
        df = vcat(df_fin, df_gow, df_lego)
        @info "Working social data..."
        Taylor.compute(df; save=true, filename=Meris.DATADIR * "macro/taylor/social.jld2")
    
        df = nothing
        df_fin = df_gow = df_lego = nothing
        GC.gc()
    
    ## BIOLOGY ##
    elseif "biology" in categories
        
        @info "Loading biology data..."
        
        #| BCI.Tree |#
        df_bci = Meris.BCITreeLoader.load()
        df_bci.class .= "eco-BCI"
        select!(df_bci, :class, :sample_id, :component_id, :counts, :nreads)
        
        #| BIOTIME |#
        df_bio = Meris.BioTIMELoader.load()
        df_bio.class .= "eco-BT" .* string.(df_bio.class)
        select!(df_bio, :class, :sample_id, :component_id, :counts, :nreads)
        
        #| GTEx |#
        df_gtex = Meris.GTExLoader.load()
        df_gtex.class .= "gen-" .* df_gtex.class
        select!(df_gtex, :class, :sample_id, :component_id, :counts, :nreads)
        
        df = vcat(df_gtex, df_bci, df_bio)
        @info "Working biology data..."
        Taylor.compute(df; save=true, filename=Meris.DATADIR * "macro/taylor/biology.jld2")
    
        df = nothing
        df_gtex = df_bci = df_bio = nothing
        GC.gc()
    end
end

### MAKE FIGURE 2 ###
function plot(;
        ext="pdf",
        big_limits=(-2, 2, -4, 4),
        small_limits=reverse([[-2,4,-3,6], [-2,5,-5,8], [-5,5,-10,10], [-6,6,-12,12]]),
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
        small_colgap=small_colgap
        )
    save(Meris.FIGDIR * "fig2.$ext", fig, pt_per_unit=1)
    return fig
end

end
