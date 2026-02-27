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

include("./../plot/colors/shadetester.jl")
using .Shades: shades

### DATA PREPARATION ###
function prepare(;
        categories=["linguistic", "microbial", "social", "biology"],
        GOFDIR=Meris.DATADIR * "goodness-of-fit/"
    )
    _compute_and_save(df, df_heaps, pareto, filename) = begin
        mkpath(dirname(filename))
        SAD.compute(
            df;
            pareto=pareto,
            nbins=30,
            heaps_df=df_heaps,
            save=true,
            filename=filename
        )
    end

    ## LINGUISTIC ##
    if "linguistic" in categories
        @info "Processing linguistic data..."
        
        #| arXiv |#
        pareto = :ParetoI
        df_full = Meris.arXivLoader.load()
        df_fit = filter_df(df_full, GOFDIR * "arxiv-candidatefits.jld2", pareto)
        df_full.class .= "arx-" .* uppercase.(df_full.class)
        df_fit.class .= "arx-" .* uppercase.(df_fit.class)
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/linguistic/arxiv.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
        GC.gc()
        
        #| Gutenberg |#
        pareto = :ParetoI
        df_full = Meris.GutenbergLoader.load()
        df_fit = filter_df(df_full, GOFDIR * "gutenberg-candidatefits.jld2", pareto)
        df_full.class .= "gutenberg-" .* uppercase.(df_full.class)
        df_fit.class .= "gutenberg-" .* uppercase.(df_fit.class)
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/linguistic/gutenberg.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
        GC.gc()
        
        #| RFCs |#
        pareto = :ParetoI
        df_full = Meris.RFCLoader.load()
        df_fit = filter_df(df_full, GOFDIR * "rfc-candidatefits.jld2", pareto)
        df_full.class .= "RFC"
        df_fit.class .= "RFC"
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/linguistic/rfc.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
        GC.gc()
    end
    
    ## MICROBIAL ##
    if "microbial" in categories
        
        @info "Processing microbial data..."
        
        #| OTU |#
        pareto = :ParetoI
        df_full = Meris.OTULoader.load()
        df_fit = filter_df(df_full, GOFDIR * "otu-candidatefits.jld2", pareto)
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/microbial/otu.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
        GC.gc()
    end

    ## SOCIAL ##
    if "social" in categories
        
        @info "Processing social data..."
        
        #| FINANCE |#
        pareto = :TemperedPareto
        df_full = Meris.FinanceLoader.load()
        df_full.class = [split(w, "-")[1] for w in df_full.class]
        df_fit = filter_df(df_full, GOFDIR * "finance-candidatefits.jld2", pareto)
        df_full.class .= "stock-" .* uppercase.(df_full.class)
        df_fit.class .= "stock-" .* uppercase.(df_fit.class)
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/social/finance.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
        GC.gc()
        
        #| Gowalla |#
        pareto = :TemperedPareto
        df_full = Meris.GowallaLoader.load()
        df_fit = filter_df(df_full, GOFDIR * "gowalla-candidatefits.jld2", pareto)
        df_full.class .= "GOWALLA"
        df_fit.class .= "GOWALLA"
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/social/gowalla.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
        GC.gc()
        
        #| LEGO |#
        pareto = :TemperedPareto
        df_full = Meris.LEGOLoader.load()
        df_fit = filter_df(df_full, GOFDIR * "lego-candidatefits.jld2", pareto)
        df_full.class .= "LEGO"
        df_fit.class .= "LEGO"
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/social/lego.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
        GC.gc()
    end

    ## BIOLOGY ##
    if "biology" in categories
        
        @info "Processing biology data..."
        
        #| BCI.Tree |#
        pareto = :ParetoI
        df_full = Meris.BCITreeLoader.load()
        df_fit = filter_df(df_full, GOFDIR * "bcitrees-candidatefits.jld2", pareto)
        df_full.class .= "eco-BCI.Trees"
        df_fit.class .= "eco-BCI.Trees"
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/biology/bci.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
        GC.gc()
        
        #| BIOTIME |#
        pareto = :ParetoI
        df_full = Meris.BioTIMELoader.load()
        df_full = df_full[df_full.class .!= "eco-BT911", :]
        df_fit = filter_df(df_full, GOFDIR * "biotime-candidatefits.jld2", pareto)
        df_full.class .= "eco-BT" .* string.(df_full.class)
        df_fit.class .= "eco-BT" .* string.(df_fit.class)
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/biology/biotime.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
        GC.gc()
        
        #| GTEx |#
        pareto = :ParetoI
        df_full = Meris.GTExLoader.load()
        df_fit = filter_df(df_full, GOFDIR * "gtex-candidatefits.jld2", pareto)
        df_full.class .= "gen-" .* string.(df_full.class)
        df_fit.class .= "gen-" .* string.(df_fit.class)
        _compute_and_save(df_fit, df_full, pareto, Meris.DATADIR * "fig3/biology/gtex.jld2")
    
        df_fit = df_full = fitdf = aicdf = nothing
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
    palette1 = shades(bases[1], 10)
    palette2 = shades(bases[2], 10)
    palette3 = shades(bases[3], 8)
    palette4 = vcat(shades(bases[4], 10)[1:6], shades(bases[5], 8))
    
    dirs = [
        Meris.DATADIR * "fig3/linguistic/",
        Meris.DATADIR * "fig3/microbial/",
        Meris.DATADIR * "fig3/social/",
        Meris.DATADIR * "fig3/biology/",
    ]
    
    # Nature max figure size: 183 mm x 170 mm (double-column).
    bigfig = Figure(
        size = (SADPlotter.NATURE_DOUBLE_WIDTH_PT, SADPlotter.NATURE_MAX_HEIGHT_PT*0.9),
        figure_padding = (4, 4, 4, 4)
    )
    
    # Fill 2×2 grid with your 4 datasets
    SADPlotter.plot!(
        bigfig[1,1]; DIR=dirs[1], palette=palette1,
        panel_id=1,
        font_scale=1.2,
        ax1limits=(nothing, nothing, 1, 3),
        ax2limits=(1e-5, 1e-1, 1e0, 1e5),
        ax3limits=(1e1, 1e8, 1e1, 1e8),
        icon_name="document.png",
        ax2_text_offset=(0.03, 3),
        ax3_text_offset=(1.0, 1.7),
        icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.1, valign=0.1)
    )
    SADPlotter.plot!(
        bigfig[1,2]; DIR=dirs[2], palette=palette2,
        panel_id=2,
        font_scale=1.2,
        ax1limits=(nothing, nothing, 1, 3),
        ax2limits=(1e-5, 1e-1, 1e0, 1e5),
        ax3limits=(1e1, 1e7, 1e1, 1e7),
        icon_name="bacteria.png",
        ax2_text_offset=(0.03, 3),
        ax3_text_offset=(1.0, 1.5),
        icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.1,  valign=0.1)
    )
    SADPlotter.plot!(
        bigfig[2,1]; DIR=dirs[3], palette=palette3,
        panel_id=3,
        font_scale=1.2,
        reverse_panel=true,
        ax1limits=(nothing, nothing, 1, 3),
        ax2limits=(1e-5, 1e-1, 1e0, 1e5),
        ax3limits=(1e1, 1e8, 1e1, 1e8),
        icon_name="socio-economic.png",
        ax2_text_offset=(0.1, 0.8),
        ax3_text_offset=(1.0, 1.7),
        icon_kw=(; width=Relative(0.77*0.25), height=Relative(0.25), halign=0.1,  valign=0.1)
    )
    SADPlotter.plot!(
        bigfig[2,2]; DIR=dirs[4], palette=palette4,
        panel_id=4,
        font_scale=1.2,
        reverse_panel=true,
        ax1limits=(nothing, nothing, 1, 3),
        ax2limits=(1e-5, 1e-1, 1e0, 1e5),
        ax3limits=(1e1, 1e7, 1e1, 1e7),
        icon_name="eco.png",
        ax2_text_offset=(0.02, 3),
        ax3_text_offset=(1.0, 1.7),
        icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.1, valign=0.1)
    )
    
    rowsize!(bigfig.layout, 1, Relative(0.43))
    rowsize!(bigfig.layout, 2, Relative(0.5))
    colsize!(bigfig.layout, 1, Relative(0.5))
    colsize!(bigfig.layout, 2, Relative(0.5))
    
    rowgap!(bigfig.layout, 0)
    colgap!(bigfig.layout, -15)
    
    save(Meris.FIGDIR * "fig3.$ext", bigfig, pt_per_unit=1)
    return bigfig
end

### HELPER ###
function filter_df(df, filename, pareto)
    fitdf = load(filename)["fitdf"]
    rename!(fitdf, :environment => :class)
    aicdf = load(filename)["aicdf"]
    rename!(aicdf, :environment => :class)
    aicdf = aicdf[(aicdf.pvalue .> 0.1) .& (aicdf.ntail .> 50), :]
    select!(aicdf, :class, :sample_id)
    fitdf = innerjoin(fitdf, aicdf, on=[:class, :sample_id])
    df = innerjoin(df, fitdf, on=[:class, :sample_id])
    select!(df, :class, :sample_id, :component_id, :counts, :nreads, pareto)
    return df
end

end # end module Figure3
