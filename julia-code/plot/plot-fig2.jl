#= Script to plot Figure 2 of main paper =#

using Pkg
Pkg.activate("./..")

using Meris
using DataFrames, DataFramesMeta, StatsBase
using CairoMakie, MakiePublication, LaTeXStrings
using JLD2

include("./../scripts/module-scripts//macropatterns/taylor.jl")
using .Taylor

include("./../plot/plot-taylor.jl")
using .TaylorPlotter

### DATA PREPARATION ###

## LINGUISTIC ##

#| arXiv |#
df_arxiv = Meris.arXivLoader.load(stopwords=true)
df_arxiv.class .= "arX-" .* df_arxiv.domain
select!(df_arxiv, [:class, :sample_id, :component_id, :counts, :nreads])

#| Gutenberg |#
df_gut = Meris.GutenbergLoader.load()
Meris.DataTools.df_filter!(df_gut; min_samples=30, min_nreads=10000, min_species=500)
df_gut.class = "guten-" .* df_gut.class
select!(df_gut, [:class, :sample_id, :component_id, :counts, :nreads])

#| RFCs |#
df_rfc = Meris.RFCLoader.load()
Meris.DataTools.df_filter!(df_rfc; min_samples=30, min_nreads=10000, min_species=500)
select!(df_rfc, [:class, :sample_id, :component_id, :counts, :nreads])

df = vcat([df_arxiv, df_gut, df_rfc]...)
tldf = Taylor.compute(df; save=true, filename=Meris.DATADIR * "macro/taylor/linguistic.jld2")

## MICROBIAL ##

#| OTU |#
df_otu = Meris.OTULoader.load()
Meris.DataTools.df_filter!(df_otu; min_samples=30, min_nreads=10000, min_species=500)
df_otu = df_otu[df_otu.class .!= "VAGINAL", :]

df = vcat([df_otu]...)
tldf = Taylor.compute(df; save=true, filename=Meris.DATADIR * "macro/taylor/microbial.jld2")

## SOCIAL ##

#| FINANCE |#
df_fin = Meris.FinanceLoader.load()
df_fin = df_fin[df_fin.class .!= "hourly",:]
df_fin.class .= "stocks-" .* df_fin.class

#| Gowalla |#
df_gow = Meris.GowallaLoader.load()
Meris.DataTools.df_filter!(df_gow; min_samples=30, min_nreads=10000, min_species=500)
df_gow.class .= "CHECK-IN"

#| LEGO |#
df_lego = Meris.LegoLoader.load(; nthemes=20)
df_lego.class .= "LEGO"
Meris.DataTools.df_filter!(df_lego; min_samples=30, min_nreads=3000, min_species=500)

df = vcat([df_fin, df_gow, df_lego]...)
tldf = Taylor.compute(df; save=true, filename=Meris.DATADIR * "macro/taylor/social.jld2")

## BIOLOGY ##

#| BCI.Tree |#
df_bci = Meris.BCITreeLoader.load(; steps=2)
df_bci.class .= "BCI.Tree"

#| GTEx |#
df_gtex = Meris.GTExLoader.load()

df = vcat([df_gtex, df_bci]...)
tldf = Taylor.compute(df; save=true, filename=Meris.DATADIR * "macro/taylor/biology.jld2")

### MAKE FIGURE 2 ###

palette1 = reverse(["#E3F2FDFF", "#BBDEFBFF", "#90CAF9FF", "#64B5F6FF", "#42A5F5FF", "#2196F3FF",
    "#1E88E5FF", "#1976D2FF", "#1565C0FF", "#0D47A1FF"])
palette2 = ["#341B0EFF", "#5B2E16FF", "#673419FF", "#79421DFF", "#A46425FF", "#B28351FF", 
    "#E2AF6DFF", "#DFB77DFF", "#FFBC38FF", "#FCCA60FF", "#F8DC8CFF"]
palette3 = reverse(["#F3CBD3FF", "#EAA9BDFF", "#DD88ACFF", "#CA699DFF", "#B14D8EFF", "#91357DFF", "#6C2167FF"])
palette4 = ["#2E7D32FF", "#388E3CFF", "#43A047FF", "#4CAF50FF", "#EF6F00FF"]

palettes = [palette1, palette2, palette3, palette4]  # 4 palettes for the 4 datasets

# Big figure
width  = 0.95 * 246
height = 3 * width / 4.67

fig = Figure(size = (2.5 * width, 1.75 * height), figure_padding = (8, 8, 8, 8))
TaylorPlotter.plot!(fig[1,1]; palettes=palettes)  # embeddable
save(Meris.FIGDIR * "fig2.pdf", fig, pt_per_unit=2)



