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
import Meris.DATADIR as DATADIR

#################
### FUNCTIONS ###
function plot_scaling(;
    αv = [0.0, 0.2, 0.5],
    DIR = DATADIR * "heap/pitmanyor/",
    BOOKDIR = DATADIR * "heap/books/",
    LEGODIR = DATADIR * "heap/lego/",
    rescale = true,
    savefig = false,
    figname = true
)
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = .85 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"\textrm{sample size}\;\log_{10}N", xlabelsize=11,
        ylabel=L"\textrm{resc. vocab. size}\;\log_{10} V(N)", ylabelsize=11,
        limits=(2,6,0,6)
    )
    #~ Define some functions for fitting
    logfun(x,p) = @. p[1] + p[2]*log(x)
    powfun(x,p) = @. p[1] + p[2]*x^p[3]

    
    # SYNTHETIC DATA / SIMULATIONS
    for i in eachindex(αv)
        αs = round(αv[i], digits=1)
        if αv[i] == 0.0
            # ~ Fit a line
            Nmin = 2
            Nmax = 6
            Nfit = exp10.(range(Nmin, Nmax, 128))
            llog = lines!(
                ax, log10.(Nfit), log10.(log.(Nfit)), color=colors[1],
                linewidth=1., label=L"\alpha=%$(αs)\;\textrm{(DP)}"
            )
        elseif αv[i] == 0.9            
            filename = "pitmanyor-vocabsize_a$(αs).jld2"
            db = JLD2.load(DIR*filename)
            V, N = db["V"], db["N"]
            # V = dropdims(mean(V, dims=2), dims=2)

            # ~ Fit a line
            Nmin = 2
            Nmax = 6
            Nfit = exp10.(range(Nmin, Nmax, 128))
            #~ Idenfity index in data that is closest to Nmin
            idx = argmin(abs.(log10.(N) .- Nmin))
            pyfit = LsqFit.curve_fit(powfun, N[idx:end], V[idx:end], [1.,1.,.5])
            if rescale
                V = @. (V - pyfit.param[1]) / pyfit.param[2]
            end
            l = lines!(
                ax, log10.(Nfit), log10.(Nfit.^αv[i]), color=colors[i+2],
                linewidth=1., label=L"\alpha=%$(αs)\;\textrm{(PY)}"
            )            
            spow = scatter!(
                ax, log10.(N[idx:end]), log10.(V[idx:end]),
                marker=:circle, markersize=5, color=:white, strokecolor=colors[i+2]
            )
        end
    end
    
    #/ DATA
    #/ Load vocabulary for books
    bookfilename = "stone-vocabsize.jld2"
    bookdb = JLD2.load(BOOKDIR*"chinese/"*bookfilename)
    Vbook, Nbook = bookdb["V"], bookdb["N"]
    Vbook = dropdims(mean(Vbook, dims=2), dims=2)
    #~ Fit a log-curve between Nmin and Nmax
    Nmin = 3.25
    Nmax = log10.(Nbook[end]) + 1
    Nfit = exp10.(range(Nmin, Nmax, 128))
    #~ Idenfity index in data that is closest to Nmin
    idx = argmin(abs.(log10.(Nbook) .- Nmin))
    
    bookfit = LsqFit.curve_fit(logfun, Nbook[idx:end], Vbook[idx:end], [1.,1.])
    if rescale
        Vbook = @. (Vbook - bookfit.param[1]) / bookfit.param[2]
    end
    #~ scatter plot raw data
    bookscat = scatter!(
        ax, log10.(Nbook), log10.(Vbook), marker=:rect, markersize=5,
        label=L"\textrm{Chinese book}", color=:white, strokecolor=colors[1]
    )    
    #~ Some text
    ct = text!(
        ax, 5.5, 0.8, color=:black,
        text = L"\propto \log N", fontsize=10, align=(:center,:top)
    )

    #/ Load vocabulary for LEGO
    legofilename = "lego-vocabsize.jld2"
    legodb = JLD2.load(LEGODIR*legofilename)
    Vlego, Nlego = legodb["V"], legodb["N"]
    Vlego = dropdims(mean(Vlego, dims=2), dims=2)
    #~ Idenfity index in data that is closest to Nmin
    Nmin = 2.0
    Nfit = exp10.(range(Nmin, Nmax, 128))
    idx = argmin(abs.(log10.(Nlego) .- Nmin))
    legofit = LsqFit.curve_fit(powfun, Nlego[idx:end], Vlego[idx:end], [1.,1.,.5])
    if rescale
        Vlego = @. (Vlego - legofit.param[1]) / legofit.param[2]
    end
    #~ scatter plot raw data
    l = lines!(
        ax, log10.(Nfit), log10.(Nfit.^legofit.param[3]),
        linewidth=1., color=colors[2], linestyle=(:dashdot,:dense)
    )
    legoscat = scatter!(
        ax, log10.(Nlego), log10.(Vlego), marker=:rect, markersize=5,
        label=L"\textrm{LEGO}", color=:white, strokecolor=colors[2]
    )
    #~ Some text
    αlego = round(legofit.param[3], digits=2)
    ct = text!(
        ax, 5.5, 2.7, color=:black,
        text = L"\propto N^{%$(αlego)}", fontsize=10, align=(:center,:bottom)
    )

    
    axislegend(
        ax,
        position=:lt, labelsize=10, patchsize=(8,20),
        margin=(4,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    
    return fig
end

end # module HeapPlotter
#/ End module
