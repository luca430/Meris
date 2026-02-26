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
    end
   
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
    end

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
    end

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
        icon_name="document.png",
        icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.1, valign=0.1)
    )
    SADPlotter.plot!(
        bigfig[1,2]; DIR=dirs[2], palette=palette2,
        ax1limits=(nothing, nothing, 1, 3),
        ax2limits=(nothing, nothing,1e0,3e4),
        ax3limits=(5, 1e7, 5, 1e5),
        icon_name="bacteria.png",
        icon_kw=(; width=Relative(0.25), height=Relative(0.25), halign=0.1,  valign=0.1)
    )
    SADPlotter.plot!(
        bigfig[2,1]; DIR=dirs[3], palette=palette3,
        reverse_panel=true,
        ax1limits=(nothing, nothing, 1, 4),
        ax3limits=(5, 5e10, 5, 5e6),
        icon_name="socio-economic.png",
        icon_kw=(; width=Relative(0.77*0.25), height=Relative(0.25), halign=0.1,  valign=0.1)
    )
    SADPlotter.plot!(
        bigfig[2,2]; DIR=dirs[4], palette=palette4,
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
