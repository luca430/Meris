#= Module to plot Figure 3 of main paper =#
module Figure3

using Meris
using DataFrames, DataFramesMeta, StatsBase
using CairoMakie, MakiePublication, LaTeXStrings
# using CSV, CodecZlib, Glob
using JLD2

include("./../scripts/module-scripts/macropatterns/SAD.jl")
using .SAD

include("./../plot/plot-SAD.jl")
using .SADPlotter

### DATA PREPARATION ###
function prepare(;
        set=["linguistic", "microbial", "social", "biology"],
        GOFDIR=Meris.DATADIR * "goodness-of-fit/"
    )
    _compute_and_save(df, pareto, filename) = begin
        mkpath(dirname(filename))
        SAD.compute(
            df;
            pareto=pareto,
            nbins=30,
            save=true,
            filename=filename
        )
    end

    ## LINGUISTIC ##
    if "linguistic" in set
        @info "Processing linguistic data..."
        
        #| arXiv |#
        pareto = :ParetoI
        df = Meris.arXivLoader.load()
        df = filter_df(df, GOFDIR * "arxiv-candidatefits.jld2", pareto)
        df.class .= "arx-" .* uppercase.(df.class)
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/linguistic/arxiv.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()
        
        #| Gutenberg |#
        pareto = :ParetoI
        df = Meris.GutenbergLoader.load()
        df = filter_df(df, GOFDIR * "gutenberg-candidatefits.jld2", pareto)
        df.class .= "gutenberg-" .* uppercase.(df.class)
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/linguistic/gutenberg.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()
        
        #| RFCs |#
        pareto = :ParetoI
        df = Meris.RFCLoader.load()
        df = filter_df(df, GOFDIR * "rfc-candidatefits.jld2", pareto)
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/linguistic/rfc.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()

    ## MICROBIAL ##
    end
    
    ## MICROBIAL ##
    if "microbial" in set
        
        @info "Processing microbial data..."
        
        #| OTU |#
        pareto = :TemperedPareto
        df = Meris.OTULoader.load()
        df = filter_df(df, GOFDIR * "otu-candidatefits.jld2", pareto)
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/microbial/otu.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()

    ## SOCIAL ##
    end

    ## SOCIAL ##
    if "social" in set
        
        @info "Processing social data..."
        
        #| FINANCE |#
        pareto = :ParetoI
        df = Meris.FinanceLoader.load()
        df = filter_df(df, GOFDIR * "finance-candidatefits.jld2", pareto)
        df.class .= "stock-" .* uppercase(df.class)
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/social/finance.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()
        
        #| Gowalla |#
        pareto = :ParetoI
        df = Meris.GowallaLoader.load()
        df = filter_df(df, GOFDIR * "gowalla-candidatefits.jld2", pareto)
        df.class .= "GOWALLA"
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/social/gowalla.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()
        
        #| LEGO |#
        pareto = :ParetoI
        df = Meris.LEGOLoader.load()
        df = filter_df(df, GOFDIR * "lego-candidatefits.jld2", pareto)
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/social/lego.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()

    ## BIOLOGY ##
    end

    ## BIOLOGY ##
    if "biology" in set
        
        @info "Processing biology data..."
        
        #| BCI.Tree |#
        pareto = :GeneralizedPareto
        df = Meris.BCITreeLoader.load()
        df = filter_df(df, GOFDIR * "bcitrees-candidatefits.jld2", pareto)
        df.class .= "eco-BCI.Trees"
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/biology/bci.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()
        
        #| BIOTIME |#
        pareto = :GeneralizedPareto
        df = Meris.BioTIMELoader.load()
        df = filter_df(df, GOFDIR * "biotime-candidatefits.jld2", pareto)
        df.class .= "eco-BT" .* string.(df.class)
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/biology/biotime.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()
        
        #| GTEx |#
        pareto = :ParetoI
        df = Meris.GTExLoader.load()
        df = filter_df(df, GOFDIR * "gtex-candidatefits.jld2", pareto)
        df.class .= "gen-" .* string.(df.class)
        _compute_and_save(df, pareto, Meris.DATADIR * "fig3/biology/gtex.jld2")
    
        df = fitdf = aicdf = nothing
        GC.gc()
    end
end

### MAKE FIGURE 2 ###
function plot(; ext="pdf")
    palette1 = reverse(["#E3F2FDFF", "#BBDEFBFF", "#90CAF9FF", "#64B5F6FF", "#42A5F5FF", "#2196F3FF",
        "#1E88E5FF", "#1976D2FF", "#1565C0FF", "#0D47A1FF"])
    palette2 = ["#341B0EFF", "#5B2E16FF", "#673419FF", "#79421DFF", "#A46425FF", "#B28351FF", 
        "#E2AF6DFF", "#DFB77DFF", "#FFBC38FF", "#FCCA60FF", "#F8DC8CFF"]
    palette3 = reverse(["#F3CBD3FF", "#EAA9BDFF", "#DD88ACFF", "#CA699DFF", "#B14D8EFF", "#91357DFF", "#6C2167FF"])
    palette4 = ["#E65100FF", "#EF6C00FF", "#F57C00FF", "#FB8C00FF", "#FF9800FF",
        "#2E7D32FF", "#388E3CFF", "#43A047FF", "#4CAF50FF"]
    
    dirs = [
        Meris.DATADIR * "fig3/microbial/",
        Meris.DATADIR * "fig3/microbial/",
        Meris.DATADIR * "fig3/biology/",
        Meris.DATADIR * "fig3/biology/",
    ]
    
    # Big figure
    width  = 0.95 * 246
    height = 3 * width / 4.67
    
    bigfig = Figure(size = (2 * 1.5*width, 2.5 * 1.5*height), figure_padding = (8, 8, 8, 8))
    
    # Fill 2×2 grid with your 4 datasets
    SADPlotter.plot!(
        bigfig[1,1]; DIR=dirs[1], palette=palette1,
        ax1limits=(nothing, nothing, 1.5, 3),
        ax2limits=(nothing, nothing, 1e1, 3e4),
        ax3limits=(5, 1e7, 5, 1e6),
        icon_name="linguistic.png",
        icon_kw=(; width=Relative(0.25), height=Relative(0.3), halign=0.05, valign=0.05)
    )
    SADPlotter.plot!(
        bigfig[1,2]; DIR=dirs[2], palette=palette2,
        ax1limits=(nothing, nothing, 1, 3),
        ax2limits=(nothing, nothing,1e0,3e4),
        ax3limits=(5, 1e7, 5, 1e5),
        icon_name="microbial.png",
        icon_kw=(; width=Relative(0.3), height=Relative(0.35), halign=0.0,  valign=0.05)
    )
    SADPlotter.plot!(
        bigfig[2,1]; DIR=dirs[3], palette=palette3,
        reverse_panel=true,
        ax1limits=(nothing, nothing, 1, 4),
        ax3limits=(5, 5e10, 5, 5e6),
        icon_name="social.png",
        icon_kw=(; width=Relative(0.4), height=Relative(0.45), halign=0.0,  valign=0.0)
    )
    SADPlotter.plot!(
        bigfig[2,2]; DIR=dirs[4], palette=palette4,
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
    
    save(Meris.FIGDIR * "fig3.$ext", bigfig, pt_per_unit=2)  # higher resolution
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
