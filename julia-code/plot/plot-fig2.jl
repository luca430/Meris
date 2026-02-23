#= Script to plot Figure 2 of main paper =#

using Meris
using DataFrames, DataFramesMeta, StatsBase
using CairoMakie, MakiePublication, LaTeXStrings
using JLD2

include("./../scripts/module-scripts/macropatterns/taylor.jl")
using .Taylor

include("./../plot/plot-taylor.jl")
using .TaylorPlotter

### DATA PREPARATION ###
function prepare()
    
    ## LINGUISTIC ##
    @info "Loading linguistic data..."
    
    #| arXiv |#
    df_arxiv = Meris.arXivLoader.load()
    df_arxiv.class .= "arx-" .* uppercase.(df_arxiv.domain)
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
    df_lego = Meris.LegoLoader.load()
    df_lego.class .= "LEGO"
    select!(df_lego, :class, :sample_id, :component_id, :counts, :nreads)
    
    df = vcat(df_fin, df_gow, df_lego)
    @info "Working social data..."
    Taylor.compute(df; save=true, filename=Meris.DATADIR * "macro/taylor/social.jld2")

    df = nothing
    df_fin = df_gow = df_lego = nothing
    GC.gc()
    
    ## BIOLOGY ##
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

### MAKE FIGURE 2 ###
fucntion plot(;ext="pdf")
    palette1 = ["#E3F2FDFF", "#BBDEFBFF", "#90CAF9FF", "#64B5F6FF", "#42A5F5FF", "#2196F3FF",
        "#1E88E5FF", "#1976D2FF", "#1565C0FF", "#0D47A1FF"]
    palette2 = reverse(["#341B0EFF", "#5B2E16FF", "#673419FF", "#79421DFF", "#A46425FF", "#B28351FF", 
        "#E2AF6DFF", "#DFB77DFF", "#FFBC38FF", "#FCCA60FF", "#F8DC8CFF"])
    palette3 = ["#F3CBD3FF", "#EAA9BDFF", "#DD88ACFF", "#CA699DFF", "#B14D8EFF", "#91357DFF", "#6C2167FF"]
    palette4 = reverse(["#E65100FF", "#EF6C00FF", "#F57C00FF", "#FB8C00FF", "#FF9800FF",
        "#2E7D32FF", "#388E3CFF", "#43A047FF", "#4CAF50FF"])
    
    palettes = [palette1, palette2, palette3, palette4]  # 4 palettes for the 4 datasets
    
    # Big figure
    width  = 0.95 * 246
    height = 3 * width / 4.67
    
    fig = Figure(size = (2.5 * width, 1.75 * height), figure_padding = (8, 8, 8, 8))
    TaylorPlotter.plot!(fig[1,1]; palettes=palettes,
        small_limits=reverse([[-2,4,-3,6], [-3,5,-5,10], [-7,6,-14,12], [-5,5,-9,10]])
        )
    save(Meris.FIGDIR * "fig2.$ext", fig, pt_per_unit=2)
end



