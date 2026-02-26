#= Module to plot Figure 3 of main paper =#
module Figure3

using Meris
using DataFrames, DataFramesMeta, StatsBase
using Colors, CairoMakie, MakiePublication
using LaTeXStrings
using JLD2

include("./../scripts/module-scripts/macropatterns/SAD.jl")
using .SAD

include("./../plot/plot-SAD.jl")
using .SADPlotter

### DATA PREPARATION ###
function prepare(
    ; set=["linguistic", "microbial", "social", "biology"],
    DIR = Meris.DATADIR * "goodness-of-fit/",
    top = 10
    )

    ## LINGUISTIC ##
    if "linguistic" in set
        @info "Loading linguistic data..."
        
        #| arXiv |#        
        minsamples    = 30      #~ min. no. of articles per "class"
        minreads      = 4_000   #~ min. no. of words per article
        mincomponents = 10_000  #~ min. no. of distinct components per "class"
        df_arxiv = Meris.arXivLoader.load(
            ; applyfilter=true, top=top,
            minsamples=minsamples, minreads=minreads, mincomponents=mincomponents
        )
        df_arxiv.class .= "arx-" .* uppercase.(df_arxiv.domain)
        select!(df_arxiv, :class, :sample_id, :component_id, :counts, :nreads)
        df_arxiv.sample_id .= string.(df_arxiv.class) .* string.(df_arxiv.sample_id)
        df_arxiv = innerjoin(
            df_arxiv, get_gof_samples(DIR *"arxiv-candidatefits.jld2"), on=[:sample_id]
        )
        
        #| Gutenberg |#
        df_gut = Meris.GutenbergLoader.load(; top=top)
        df_gut.class = "guten-" .* uppercase.(df_gut.class)
        select!(df_gut, :class, :sample_id, :component_id, :counts, :nreads)
        df_gut.sample_id .= string.(df_gut.class) .* string.(df_gut.sample_id)
        df_gut = innerjoin(
            df_gut, get_gof_samples(DIR * "gutenberg-candidatefits.jld2"), on=[:sample_id]
        )
        
        #| RFCs |#
        df_rfc = Meris.RFCLoader.load(; top=top)
        df_rfc.class .= uppercase.(df_rfc.class)
        select!(df_rfc, :class, :sample_id, :component_id, :counts, :nreads)
        df_rfc.sample_id .= string.(df_rfc.class) .* string.(df_rfc.sample_id)
        df_rfc = innerjoin(
            df_rfc, get_gof_samples(DIR * "rfc-candidatefits.jld2"), on=[:sample_id]
        )
        
        df = vcat(df_arxiv, df_gut, df_rfc)
        @info "Preparing test/linguistic data..."
        SAD.compute(
                df;
                xmins=10 .^ collect(-5.5:0.01:-3.5),
                pareto_type="I",
                nbins=25,
                filter=true,
                save=true,
                filename=Meris.DATADIR * "macro/sad/linguistic.jld2"
            )
    
        df = nothing
        df_arxiv = df_gut = df_rfc = nothing
        GC.gc()
    end
   
    ## MICROBIAL ##
    if "microbial" in set        
        @info "Loading microbial data..."
        
        #| OTU |#
        df_otu = Meris.OTULoader.load()
        select!(df_otu, :class, :sample_id, :component_id, :counts, :nreads)
        df_otu.sample_id .= string.(df_otu.class) .* string.(df_otu.sample_id)
        df_otu = innerjoin(df_otu, get_gof_samples(Meris.DATADIR * "gof/otu-candidatefits.jld2"), on=[:sample_id])
        
        df = vcat(df_otu)
        @info "Preparing microbial data..."
        SAD.compute(
                df;
                xmins=10 .^ collect(-5.0:0.01:-3.5),
                pareto_type="I",
                nbins=25,
                filter=true,
                save=true,
                filename=Meris.DATADIR * "macro/sad/microbial.jld2"
            )
    
        df = nothing
        df_otu = nothing
        GC.gc()
    end

    ## SOCIAL ##
    if "social" in set        
        @info "Loading social data..."
        
        #| FINANCE |#
        df_fin = Meris.FinanceLoader.load()
        df_fin = df_fin[endswith.(df_fin.class, "-daily"), :]
        df_fin.class = replace.(df_fin.class, "-daily" => "")
        df_fin.class .= "stock-" .* uppercase.(df_fin.class)
        select!(df_fin, :class, :sample_id, :component_id, :counts, :nreads)
        df_fin.sample_id .= string.(df_fin.class) .* string.(df_fin.sample_id)
        df_fin = innerjoin(df_fin, get_gof_samples(Meris.DATADIR * "gof/finance-candidatefits.jld2"), on=[:sample_id])
        
        #| Gowalla |#
        df_gow = Meris.GowallaLoader.load()
        df_gow.class .= "CHECK-IN"
        select!(df_gow, :class, :sample_id, :component_id, :counts, :nreads)
        df_gow.sample_id .= string.(df_gow.class) .* string.(df_gow.sample_id)
        df_gow = innerjoin(df_gow, get_gof_samples(Meris.DATADIR * "gof/gowalla-candidatefits.jld2"), on=[:sample_id])
        
        #| LEGO |#
        df_lego = Meris.LegoLoader.load()
        df_lego.class .= "LEGO"
        select!(df_lego, :class, :sample_id, :component_id, :counts, :nreads)
        df_lego.sample_id .= string.(df_lego.class) .* string.(df_lego.sample_id)
        
        df = vcat(df_fin, df_gow, df_lego)
        @info "Preparing social data..."
        SAD.compute(
                df;
                xmins=10 .^ collect(-5.5:0.01:-2.0),
                pareto_type="I",
                nbins=25,
                filter=true,
                save=true,
                filename=Meris.DATADIR * "macro/sad/social.jld2"
            )
        
    
        df = nothing
        df_fin = df_gow = df_lego = nothing
        GC.gc()
    end

    ## BIOLOGY ##
    if "biology" in set        
        @info "Loading biology data..."
        
        #| BCI.Tree |#
        df_bci = Meris.BCITreeLoader.load()
        df_bci.class .= "eco-BCI"
        select!(df_bci, :class, :sample_id, :component_id, :counts, :nreads)
        df_bci.sample_id .= string.(df_bci.class) .* string.(df_bci.sample_id)
        df_bci = innerjoin(df_bci, get_gof_samples(Meris.DATADIR * "gof/bcitrees-candidatefits.jld2"), on=[:sample_id])
        
        #| BIOTIME |#
        df_bio = Meris.BioTIMELoader.load()
        df_bio.class .= "eco-BT" .* string.(df_bio.class)
        select!(df_bio, :class, :sample_id, :component_id, :counts, :nreads)
        df_bio.sample_id .= string.(df_bio.class) .* string.(df_bio.sample_id)
        df_bio = innerjoin(df_bio, get_gof_samples(Meris.DATADIR * "gof/biotime-candidatefits.jld2"), on=[:sample_id])
        
        #| GTEx |#
        df_gtex = Meris.GTExLoader.load()
        df_gtex.class .= "gen-" .* df_gtex.class
        select!(df_gtex, :class, :sample_id, :component_id, :counts, :nreads)
        df_gtex.sample_id .= string.(df_gtex.class) .* string.(df_gtex.sample_id)
        df_gtex = innerjoin(df_gtex, get_gof_samples(Meris.DATADIR * "gof/gtex-candidatefits.jld2"), on=[:sample_id])
        
        df = vcat(df_gtex, df_bci, df_bio)
        @info "Preparing biology data..."
        SAD.compute(
                df;
                xmins=10 .^ collect(-5.0:0.01:-2.0),
                pareto_type="I",
                nbins=25,
                filter=true,
                save=true,
                filename=Meris.DATADIR * "macro/sad/biology.jld2"
            )
    
        df = nothing
        df_gtex = df_bci = df_bio = nothing
        GC.gc()
    end
end

### MAKE FIGURE 2 ###
function plot(
    ; ext="pdf",
    bases = [
        colorant"#1f77b4",  # blue
        colorant"#ff7f0e",  # orange
        colorant"#9467bd",  # purple
        colorant"#2ca02c",  # green
        colorant"#d62728"   # red
    ]
    )
    palette1 = shades(bases[1], 7)
    palette2 = shades(bases[2], 10)
    palette3 = shades(bases[3], 5)
    palette4 = vcat(shades(bases[4], 5), shades(bases[5], 4))
    
    zipfdirs = [
        Meris.DATADIR * "macro/sad/linguistic.jld2",
        Meris.DATADIR * "macro/sad/microbial.jld2",
        Meris.DATADIR * "macro/sad/social.jld2",
        Meris.DATADIR * "macro/sad/microbial.jld2",
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
        icon_name="document.png",
        icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.1, valign=0.1)
    )
    SADPlotter.plot!(
        bigfig[1,2]; ZIPFDIR=zipfdirs[2], palette=palette2,
        ax1limits=(nothing, nothing, 1, 3),
        ax2limits=(nothing, nothing,1e0,3e4),
        ax3limits=(5, 1e7, 5, 1e5),
        icon_name="bacteria.png",
        icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.1,  valign=0.1)
    )
    SADPlotter.plot!(
        bigfig[2,1]; ZIPFDIR=zipfdirs[3], palette=palette3,
        reverse_panel=true,
        ax1limits=(nothing, nothing, 1, 4),
        ax3limits=(5, 5e10, 5, 5e6),
        icon_name="socio-economic.png",
        icon_kw=(; width=Relative(0.77*0.25), height=Relative(0.25), halign=0.1,  valign=0.1)
    )
    SADPlotter.plot!(
        bigfig[2,2]; ZIPFDIR=zipfdirs[4], palette=palette4,
        reverse_panel=true,
        ax1limits=(nothing, nothing, 1, 3),
        ax2limits=(nothing, nothing,1e0,1e5),
        ax3limits=(5e0, 5e8, 5, 5e5),
        icon_name="eco.png",
        icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.1, valign=0.1)
    )
    
    rowsize!(bigfig.layout, 1, Relative(0.43))
    rowsize!(bigfig.layout, 2, Relative(0.5))
    colsize!(bigfig.layout, 1, Relative(0.5))
    colsize!(bigfig.layout, 2, Relative(0.5))
    
    rowgap!(bigfig.layout, 15)
    colgap!(bigfig.layout, -15)
    
    save(Meris.FIGDIR * "fig3.$ext", bigfig, pt_per_unit=2)  # higher resolution
    # return bigfig
end

########################
### HELPER FUNCTIONS ###
function shades(base::Colorant, n)
    hsl_base = convert(HSL, base)
    h, s = hsl_base.h, hsl_base.s
    _colors = [HSL(h, s, l) for l in range(0.3, 0.7, length=n)]
    return _colors
end


function get_gof_samples(file)
    df = load(file)["aicdf"]
    select!(df, :environment, :sample_id)
    df.sample_id .= string.(df.environment) .* string.(df.sample_id)
    return df
end

end # end module Figure3


