#= Simple module for plotting and visualizing distributions =#
#/ Start module
module DistPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings

using Distributions
using SpecialFunctions
using LsqFit

import Meris.TweedieDistribution as Tw

#################
### FUNCTIONS ###
function plot_pareto(;
    α = 1.0,
    θ = 1.0,
    savefig = false,
    figname = nothing
)
    __theme = MakiePublication.theme_acs(; ishollowmarkers=[true,true])
    set_theme!(__theme)
    
	  width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"\log x", xlabelsize=11,
        ylabel=L"\log p(x)", ylabelsize=11,
        limits=(0,3,-6,0)
    )

    xplot = exp.(range(0,4,128))    
    
    #~ Plot Pareto
    pplot = pareto.(xplot; α=α, θ=θ)
    l = lines!(
        ax, log.(xplot), pplot, label=L"\textrm{Pareto(x; \alpha,\theta)}",
        color=:black, linewidth=.8
    )

    #~ Plot Lomax
    splot = lomax.(xplot; α=α, θ=θ)
    s = lines!(ax, log.(xplot), splot, label=L"\textrm{Lomax(x; \alpha,\theta)}", linewidth=1.)
    #~ Plot re-shifted Lomax
    rplot = lomax.(xplot .- θ; α=α, θ=θ)
    r = lines!(
        ax, log.(xplot), rplot, linewidth=1., linestyle=(:dash, :dense),
        label=L"\textrm{Lomax(x - \theta; \alpha,\theta)}", 
    )

    axislegend(
        ax,
        position=:rt, labelsize=10, patchsize=(8,20),
        margin=(8,0,0,0), patchlabelgap=2, padding=(0,0,0,0)
    )
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

function plot_TruncatedPareto4(;
    α = 1.5,
    β = 1.0,
    γ = 1.0,
    ε = 1.0,
    savefig = false,
    figname = nothing
)
    __theme = MakiePublication.theme_acs(; ishollowmarkers=[true,true])
    set_theme!(__theme)
    
	  width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"\log x", xlabelsize=11,
        ylabel=L"\log p(x)", ylabelsize=11,
    )

    xplot = exp.(range(-4,4,128))
    
    #~ Plot Pareto
    θ = [α, β, γ, ε]
    pplot = TruncatedPareto4.(xplot, Ref(θ))
    l = lines!(
        ax, log.(xplot), pplot, label=L"\textrm{Pareto(x; \alpha,\theta)}",
        color=:black, linewidth=.8
    )
    
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

"""Compare Gamma and tempered Pareto distribution"""
function plot_compare(; fitgamma=false)
    __theme = MakiePublication.theme_acs(; ishollowmarkers=[true,true])
    set_theme!(__theme)
    
    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"x", xlabelsize=11,
        ylabel=L"p(x)", ylabelsize=11,
        limits=(0,5,-30,2)
    )

    Pf(x,α,ε) = α*ε^α * x^(-α-1)
	  TePf(x,α,β,ε) = (x >= ε) ? ε^α * exp(β*ε) * x^(-α-1) * exp(-β*x) * (α + β*x) : 0.0
    TrPf(x,α,ε,κ) = (ε <= x <= κ) ? α*ε^α * x^(-α-1) / (1 - (ε/κ)^α) : 0.0
    Gf(x,α,θ) = (x > 0.0) ? x^(α-1) * exp(-x/θ) / SpecialFunctions.gamma(α) / θ^α : 0.0

    α = 0.5
    β = 1 / 1e3
    ε = 1e0
    
    x = exp10.(range(log10(ε),log10(1e5),128))
    py = Pf.(x, Ref(α), Ref(ε))
    tepy = TePf.(x, Ref(α), Ref(β), Ref(ε))
    trpy = TrPf.(x, Ref(α), Ref(ε), Ref(1e4))

    #~ Sample and fit Gamma, as normalization constants are way off
    gy = Gf.(x, Ref(1-α), Ref(1/β))

    #~ Curve fit a Gamma
    lb = [0.0, 0.0]
    ub = [Inf, Inf]

    lines!(ax, log10.(x), log.(py), label=L"\textrm{Pareto}")
    lines!(ax, log10.(x), log.(tepy), label=L"\textrm{tempPareto}")
    lines!(ax, log10.(x), log.(trpy), label=L"\textrm{truncPareto}")
    lines!(ax, log10.(x), log.(gy), label=L"\textrm{Gamma}")
    if fitgamma
        fGamma(x,p) = Gf.(x, Ref(p[1]), Ref(p[2]))
        gfit = LsqFit.curve_fit(fGamma, x, tpy, [1-α, 1/β], lower=lb, upper=ub)
        gfity = fGamma(x, gfit.param)
        lines!(ax, log10.(x), log.(gfity), label=L"\textrm{fitted Gamma}")
    end

    axislegend(ax, position=:lb, patchsize=(8,8), rowgap=0., labelsize=10, padding=0)
    return fig
end

"Plot Tweedie distribution"
function plot_tweedie(;
    μ::Float64 = 2.5,
    ϕ::Float64 = .75,
    b::Array{Float64} = [1.9, 2.6]
)
    __theme = MakiePublication.theme_acs(; ishollowmarkers=[true,true])
    set_theme!(__theme)
    
    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(1.8 .* width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"x", xlabelsize=14,
        ylabel=L"\log_{10}\;p(x)", ylabelsize=14,
        limits=(0,20,-6,0)        
    )
    logax = Axis(
        fig[1,2],
        xlabel=L"\log_{10}\;x", xlabelsize=14,
        ylabel=L"\log_{10}\;p(x)", ylabelsize=14,
        limits=(-2,2,-5,0)
    )
    x = collect(range(0.0,20.0,128))
    expx = collect(exp10.(range(-1.5,2,128)))
    
    #~ Gamma
    gm = Distributions.Gamma(1/ϕ, ϕ*μ)
    ygm = Distributions.pdf.(gm, x)
    logygm = expx .* Distributions.pdf.(gm, expx)
    lines!(
        ax, x, log10.(ygm), label=L"\textrm{gamma}",
        linewidth=.8, linestyle=(:dash,:dense), color=:black
    )
    lines!(
        logax, log10.(expx), log10.(logygm),
        linewidth=.8, linestyle=(:dash,:dense), color=:black
    )
    #~ Inv Gaussian
    invg = Distributions.InverseGaussian(μ, 1/ϕ)
    yig = Distributions.pdf.(invg, x)
    logyig = expx .* Distributions.pdf.(invg, expx)
    lines!(
        ax, x, log10.(yig), label=L"\textrm{inverse Gaussian}",
        linewidth=.8, linestyle=(:dash,:dense), color=:gray
    )
    lines!(
        logax, log10.(expx), log10.(logyig),
        linewidth=.8, linestyle=(:dash,:dense), color=:gray
    )
    #~ Poisson
    # ps = Distributions.Poisson(μ)
    # xps = 0:15
    # yps = Distributions.pdf.(ps, xps)
    # scatterlines!(
    #     ax, xps, log10.(yps), label=L"\textrm{Poisson}", color=:black,
    #     linestyle=(:dash,:dense), linewidth=.7, markersize=5
    # )

    for i in eachindex(b)
        tw = Tw.Tweedie(b[i], μ, ϕ)        
        ytw = Tw.pdf.(tw, x)
        logytw = expx .* Tw.pdf.(tw, expx)
        lines!(ax, x, log10.(ytw), label=L"\textrm{Tweedie}(b=%$(b[i]))")
        lines!(logax, log10.(expx), log10.(logytw))
    end

    axislegend(ax, position=:lb, patchsize=(12,5), rowgap=0., labelsize=12, padding=0)
    return fig
end

########################
### HELPER FUNCTIONS ###
function pareto(x; α=2.0, θ=1.0)
    p = (x < θ) ? NaN : log.(α*θ^α / x^(1+α))
    return p
end

function lomax(x; α=2.0, θ=1.0)
	  p = (x < 0) ? NaN : log.(α*θ^α / (x + θ)^(1+α))
    return p
end

function TruncatedPareto4(x, θ)
    α, β, γ, ε = θ
    p = (x < 0) ? 0.0 : x^(-α-1)*(x+ε)^(β-α-1)*exp(-γ*x/ε)
    return log(p)
end

end # module DistPlotter
#/ End module
