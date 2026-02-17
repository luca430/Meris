#= Script to plot Figure 3 of main paper =#

using Pkg
Pkg.activate("./..")

using Meris
using DataFrames, DataFramesMeta, StatsBase
using CairoMakie, MakiePublication, LaTeXStrings
using CSV, CodecZlib, Glob

include("./../scripts/module-scripts//macropatterns/SAD.jl")
using .SAD

include("./../plot/plot-SAD.jl")
using .SADPlotter

### DATA PREPARATION ###

## LINGUISTIC ##

#| arXiv |#
df_arxiv = Meris.arXivLoader.load(stopwords=true)
df_arxiv.class .= "arx-" .* df_arxiv.domain
Meris.DataTools.df_filter!(df_arxiv; min_samples=30, min_nreads=10000, min_species=100)
select!(df_arxiv, [:class, :sample_id, :component_id, :counts, :nreads])

#| Gutenberg |#
df_gut = Meris.GutenbergLoader.load()
Meris.DataTools.df_filter!(df_gut; min_samples=30, min_nreads=Int(1e5), min_species=100)
df_gut.class = "guten-" .* df_gut.class
select!(df_gut, [:class, :sample_id, :component_id, :counts, :nreads])

#| RFCs |#
df_rfc = Meris.RFCLoader.load()
Meris.DataTools.df_filter!(df_rfc; min_samples=30, min_nreads=10000, min_species=100)
select!(df_rfc, [:class, :sample_id, :component_id, :counts, :nreads])

df = vcat([df_arxiv, df_gut, df_rfc]...)
out = SAD.compute(
        df;
        xmins=10 .^ collect(-5.5:0.01:-3.5),
        pareto_type="I",
        nbins=25,
        filter=true,
        save=true,
        filename=Meris.DATADIR * "macro/sad/linguistic.jld2"
    );

## MICROBIAL ##

#| OTU |#
df_otu = Meris.OTULoader.load()
Meris.DataTools.df_filter!(df_otu; min_samples=30, min_nreads=10000, min_species=100)
df_otu = df_otu[df_otu.class .!= "VAGINAL", :]

df = vcat([df_otu]...)
out = SAD.compute(
        df;
        xmins=10 .^ collect(-5.0:0.01:-3.5),
        pareto_type="I",
        nbins=25,
        filter=true,
        save=true,
        filename=Meris.DATADIR * "macro/sad/microbial.jld2"
    );

## SOCIAL ##

#| FINANCE |#
df_fin = Meris.FinanceLoader.load()
df_fin = df_fin[endswith.(df_fin.class, "-daily"), :]
df_fin.class = replace.(df_fin.class, "-daily" => "")
df_fin.class .= "stock-" .* df_fin.class
Meris.DataTools.df_filter!(df_fin; min_samples=30, min_nreads=10000, min_species=100)

#| Gowalla |#
df_gow = Meris.GowallaLoader.load()
Meris.DataTools.df_filter!(df_gow; min_samples=30, min_nreads=10000, min_species=100)
df_gow.class .= "CHECK-IN"

#| LEGO |#
df_lego = Meris.LegoLoader.load(; nthemes=100)
df_lego.class .= "LEGO"
Meris.DataTools.df_filter!(df_lego; min_samples=30, min_nreads=3000, min_species=100)

df = vcat([df_fin, df_gow, df_lego]...)
out = SAD.compute(
        df;
        xmins=10 .^ collect(-5.5:0.01:-2.0),
        pareto_type="I",
        nbins=25,
        filter=true,
        save=true,
        filename=Meris.DATADIR * "macro/sad/social.jld2"
    );

## BIOLOGY ##

#| BCI.Tree |#
df_bci = Meris.BCITreeLoader.load(; steps=2)
df_bci.class .= "eco-BCI"
Meris.DataTools.df_filter!(df_bci; min_samples=30, min_nreads=1000, min_species=50)

#| BIOTIME |#
df_bio = Meris.BioTIMELoader.load()
df_bio.class .= "eco-BT" .* string.(df_bio.class)
df_bio = df_bio[df_bio.class .!= "eco-BT634", :]
Meris.DataTools.df_filter!(df_bio; min_samples=30, min_nreads=5000, min_species=100)

#| GTEx |#
df_gtex = Meris.GTExLoader.load()
df_gtex.class .= "gen-" .* df_gtex.class
Meris.DataTools.df_filter!(df_gtex; min_samples=30, min_nreads=Int(1e8), min_species=100)

df = vcat([df_gtex, df_bci]...)
out = SAD.compute(
        df;
        xmins=10 .^ collect(-5.0:0.01:-2.5),
        pareto_type="I",
        nbins=25,
        filter=true,
        save=true,
        filename=Meris.DATADIR * "macro/sad/biology.jld2"
    );

### MAKE FIGURE 2 ###
palette1 = reverse(["#E3F2FDFF", "#BBDEFBFF", "#90CAF9FF", "#64B5F6FF", "#42A5F5FF", "#2196F3FF",
    "#1E88E5FF", "#1976D2FF", "#1565C0FF", "#0D47A1FF"])
palette2 = ["#341B0EFF", "#5B2E16FF", "#673419FF", "#79421DFF", "#A46425FF", "#B28351FF", 
    "#E2AF6DFF", "#DFB77DFF", "#FFBC38FF", "#FCCA60FF", "#F8DC8CFF"]
palette3 = reverse(["#F3CBD3FF", "#EAA9BDFF", "#DD88ACFF", "#CA699DFF", "#B14D8EFF", "#91357DFF", "#6C2167FF"])
palette4 = ["#E65100FF", "#EF6C00FF", "#F57C00FF", "#FB8C00FF", "#FF9800FF",
    "#2E7D32FF", "#388E3CFF", "#43A047FF", "#4CAF50FF"]

# Example data sources (adjust paths)
zipfdirs = [
    Meris.DATADIR * "macro/sad/linguistic.jld2",
    Meris.DATADIR * "macro/sad/microbial.jld2",
    Meris.DATADIR * "macro/sad/social.jld2",
    Meris.DATADIR * "macro/sad/biology.jld2",
]

# Big figure
width  = 0.95 * 246
height = 3 * width / 4.67

bigfig = Figure(size = (2 * 1.5*width, 2.5 * 1.5*height), figure_padding = (8, 8, 8, 8))

# Fill 2×2 grid with your 4 datasets
SADPlotter.plot!(
    bigfig[1,1]; ZIPFDIR=zipfdirs[1], palette=palette1,
    ax1limits=(nothing, nothing, 1.5, 3),
    ax2limits=(nothing, nothing, 1e1, 3e4),
    ax3limits=(5, 1e7, 5, 1e6),
    icon_name="linguistic.png",
    icon_kw=(; width=Relative(0.25), height=Relative(0.3), halign=0.05, valign=0.05)
)
SADPlotter.plot!(
    bigfig[1,2]; ZIPFDIR=zipfdirs[2], palette=palette2,
    ax1limits=(nothing, nothing, 1, 3),
    ax2limits=(nothing, nothing,1e0,3e4),
    ax3limits=(5, 1e7, 5, 1e5),
    icon_name="microbial.png",
    icon_kw=(; width=Relative(0.3), height=Relative(0.35), halign=0.0,  valign=0.05)
)
SADPlotter.plot!(
    bigfig[2,1]; ZIPFDIR=zipfdirs[3], palette=palette3,
    reverse_panel=true,
    ax1limits=(nothing, nothing, 1, 4),
    ax3limits=(5, 5e10, 5, 5e6),
    icon_name="social.png",
    icon_kw=(; width=Relative(0.4), height=Relative(0.45), halign=0.0,  valign=0.0)
)
SADPlotter.plot!(
    bigfig[2,2]; ZIPFDIR=zipfdirs[4], palette=palette4,
    reverse_panel=true,
    ax1limits=(nothing, nothing, 1, 3),
    ax2limits=(nothing, nothing,1e0,1e5),
    ax3limits=(5e0, 5e8, 5, 5e5),
    icon_name="biology.png",
    icon_kw=(; width=Relative(0.45), height=Relative(0.5), halign=0.03, valign=-0.1)
)

rowsize!(bigfig.layout, 1, Relative(0.43))
rowsize!(bigfig.layout, 2, Relative(0.5))
colsize!(bigfig.layout, 1, Relative(0.5))
colsize!(bigfig.layout, 2, Relative(0.5))

rowgap!(bigfig.layout, 15)
colgap!(bigfig.layout, -15)

save(Meris.FIGDIR * "fig3.pdf", bigfig, pt_per_unit=2)  # higher resolution




