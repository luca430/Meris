#= Module to plot Heap's law for some processes and/or data =#
#/ Start module
module HeapPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings
using LsqFit
using JLD2
using StatsBase

#/ Modules
import Moira.DATADIR as DATADIR

#################
### FUNCTIONS ###
function plot_scaling(;
    αv = [0.0, 0.2, 0.5, 0.9],
    DIR = DATADIR * "heap/pitmanyor/",
    BOOKDIR = DATADIR * "heap/books/",
    savefig = false,
    figname = true
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    width = .85 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(width,height), figure_padding=(2,6,1,6))
    ax = Axis(
        fig[1,1],
        xlabel=L"\textrm{sample size}\;\log_{10}N", xlabelsize=11,
        ylabel=L"\textrm{vocab. size}\;\log_{10} V(N)", ylabelsize=11,
        limits=(2,6,1,5)
    )

    #/ Do the same for a Chinese book
    cfilename = "stone-vocabsize.jld2"
    cdb = JLD2.load(BOOKDIR*"chinese/"*cfilename)
    cV, cN = cdb["V"], cdb["N"]
    cV = dropdims(mean(cV, dims=2), dims=2)

    
    
    Nmin = 3.25
    Nmax = log10.(cN[end]) + 1
    Nfit = exp10.(range(Nmin, Nmax, 128))
    cs = scatter!(ax, log10.(cN), log10.(cV), marker=:rect, markersize=5)
    clogfun(x,p) = @. p[1] + p[2]*log(x)
    cfit = LsqFit.curve_fit(clogfun, cN[5:end], cV[5:end], [0.,1.,1.])
    l = lines!(
        ax, log10.(Nfit), log10.(clogfun(Nfit, cfit.param)),
        linewidth=1., label=L"\textrm{Story of the stone}", linestyle=(:dash,:dense)
    )
    tidx = 82
    ct = text!(
        ax, log10(Nfit[tidx]), log10(clogfun(Nfit[tidx], cfit.param))*1.05,
        text = L"\propto \log N", fontsize=10, align=(:right,:bottom)
    )

    for i in eachindex(αv)
        αs = round(αv[i], digits=1)
        filename = "pitmanyor-vocabsize_a$(αs).jld2"
        db = JLD2.load(DIR*filename)
        V, N = db["V"], db["N"]
        V = dropdims(mean(V, dims=2), dims=2)

        #~ Fit a line
        Nmin = 3
        Nmax = 6
        Nfit = exp10.(range(Nmin, Nmax, 128))

        fn(x, p) = (αv[i] == 0.0) ? (@. p[1] + p[2]*log(x)) : (@. p[1] + p[2] * x^αv[i])
        
        # if αv[i] == 0.0
        #     logfun(x,p) = @. p[1] + p[2]*log(x)
        #     #~ Idenfity index in data that is closest to Nmin
        #     idx = argmin(abs.(log10.(N) .- Nmin))
        #     fit = LsqFit.curve_fit(logfun, N[idx:end], V[idx:end], [0.,1.])
        #     l = lines!(
        #         ax, log10.(Nfit), log10.(logfun(Nfit, fit.param)),
        #         linewidth=1., linestyle=:dash
        #     )
        # else
        #     powlaw(x,p) = @. p[1] + p[2] * x^αv[i]
        #     #~ Idenfity index in data that is closest to Nmin
        idx = argmin(abs.(log10.(N) .- Nmin))
        fit = LsqFit.curve_fit(fn, N[idx:end], V[idx:end], [0.,1.])
        
        l = lines!(
            ax, log10.(Nfit), log10.(fn(Nfit, fit.param)),
            linewidth=1., linestyle=:dash
        )
        # end

        #~ Scatter data
        s = scatter!(ax, log10.(N), log10.(V), label=L"\alpha=%$(αs)", markersize=4.)
    end    

    
    axislegend(
        ax,
        position=:lt, labelsize=10, patchsize=(6,20),
        margin=(4,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    
    return fig
end

end # module HeapPlotter
#/ End module
