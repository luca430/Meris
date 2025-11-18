#= Simple module for plotting and visualizing distributions =#
#/ Start module
module DistPlotter

#/ Packages
using CairoMakie
using MakiePublication
using LaTeXStrings

using Distributions
using SpecialFunctions

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
function plotcompare()
    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(width,height), figure_padding=(2,4,2,14))
    ax = Axis(
        fig[1,1],
        xlabel=L"x", xlabelsize=11,
        ylabel=L"p(x)", ylabelsize=11,
        limits=(0,5,-30,2)
    )
    
	  TPf(x,α,β,ε) = ε^α * exp(β*ε) * x^(-α-1) * exp(-β*x) * (α + β*x)
    Gf(x,α,θ) = x^(α-1) * exp(-x/θ) / SpecialFunctions.gamma(α) / θ^α

    α = 1.0
    β = 1 / 1e4
    ε = 1e0
    
    x = exp10.(range(log10(ε),log10(1e5),64))
    tpy = TPf.(x, Ref(α), Ref(β), Ref(ε))

    #~ Sample and fit Gamma, as normalization constants are way off
    gy = Gf.(x, Ref(α), Ref(1/β))
    Δy = abs.(tpy .- gy)

    lines!(ax, log10.(x), log.(tpy), label=L"\textrm{tPareto}")
    lines!(ax, log10.(x), log.(gy), label=L"\textrm{Gamma}")
    lines!(ax, log10.(x), log.(Δy), label=L"\Delta", linewidth=.8, linestyle=:dash, color=:black)

    axislegend(ax, position=:lb, patchsize=(8,8), rowgap=0., labelsize=10, padding=0)
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
