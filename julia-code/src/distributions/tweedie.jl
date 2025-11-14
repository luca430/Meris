#= Tweedie distributions =#
using Distributions
using Random
using SpecialFunctions

"""
    Tweedie(p,μ,ϕ)

The *Tweedie distribution* with power exponent `p`, mean `μ` and dispersion parameter `ϕ`
has a probability density function

```math
f_p(y;\\mu,\\theta) =
\\begin{cases}
a(y,\\theta) \\exp \\left[ \\frac{1}{\\phi} \\left( y\\theta - \\kappa(\\theta) \\right) \\right]
\\end{cases}
```

```julia
Tweedie(p)      # Tweedie distribution with unit mean and dispersion, i.e. Tweedie(p,1,1)
Tweedie(p,μ,θ)  # Tweedie distribution with exponent p, mean μ and dispersion parameter ϕ

params(d)

External links

* [Tweedie distribution on Wikipedia](https://en.wikipedia.org/wiki/Tweedie_distribution)
* [Functions of a statistical variate with given means, with special reference to Laplacian distributions](https://doi.org/10.1017/S0305004100023185)
* [Series evaluation of Tweedie exponential dispersion model densities](https://doi.org/10.1007/s11222-005-4070-y)
"""
struct Tweedie{T<:Real} <: ContinuousUnivariateDistribution
    p::T
    μ::T
    ϕ::T

    Tweedie{T}(p::T, μ::T, ϕ::T) where {T<:Real} = new{T}(p, μ, ϕ)
end

function Tweedie(p::T, μ::T, ϕ::T; check_args::Bool=true) where {T<:Real}
    #~ Special cases
    (p == 0) && (return Normal(μ, sqrt(ϕ)))
    (p == 1) && (return Poisson(μ))
    (p == 2) && (return Gamma(1/ϕ,ϕ*μ))
    # (p == 3) && (return InverseGaussian(μ,1/ϕ))
	  Distributions.@check_args Tweedie (p, !(zero(p) < p < one(p))) (μ, μ >= zero(μ)) (ϕ, ϕ > zero(ϕ))
    return Tweedie{T}(p, μ, ϕ)
end

### Constructors
Tweedie(p::Real; check_args::Bool=true) = Tweedie(p, one(p), one(p); check_args=check_args)
Tweedie(p::Real, μ::Real, ϕ::Real; check_args::Bool=true) = Tweedie(promote(p, μ, ϕ)...; check_args=check_args)

### Parameters
params(d::Tweedie) = (d.p, d.μ, d.ϕ)

### Statistics
mean(d::Tweedie) = d.μ
var(d::Tweedie) = d.μ^d.p
shape(d::Tweedie) = (2 - d.p) / (d.p - 1)
θ(d::Tweedie) = isone(d.p) ? log(d.μ) : d.μ^(1 - d.p) / (1 - d.p)
κ(d::Tweedie) = isone(d.p - 1) ? log(d.μ) : d.μ^(2 - d.p) / (2 - d.p)

### Evaluation
"""
    pdf(d::Tweedie, x::Real)

The density function of a Tweedie distribution generally does not have a closed form (depending
on the value of `p`), and hence the density at `x` needs to be approximated. Here, the series
expansion of Dunn & Smyth (2005) is used. In general, the series expansion is good for values
`x` not too small (a rough estimate is x>0.1) and is excellent for x≥0.5. Simultaneously, the
series expansion, and also the maximum likelihood estimate, becomes worse as `p` approaches the
Poisson limit at p=1.
"""

function pdf(d::Tweedie, x::Real; ε::Float64 = 1e-16)
    (1 < d.p < 2) && (return _pdfpoissongamma(d::Tweedie, x::Real; ε=ε))
    (d.p > 2) && (return _pdfpoissonstable(d::Tweedie, x::Real; ε=ε))
end

"""
Approximate the infinite series W(x;ϕ,p) = ∑ₖ Wₖ using methods from Dunn & Smyth (2005)
"""
function _pdfpoissongamma(d::Tweedie, x::Real; ε::Float64=1e-16)
	  if iszero(x)
        return exp(-(d.μ^(2 - d.p)) / (d.ϕ * (2 - d.p)))
    end
    #~ Compute estimates of kmax and Wmax
    kmax = round(Int, x^(2 - d.p) / (d.ϕ * (2 - d.p)))
    α = -shape(d)
    logWmax = kmax*(α-1) - log(2π) - log(kmax) - log(-α) / 2
    #~ Walk left and right until Wₖ is less than log(Wmax) - log(ε)
end

function _pdfpoissonstable(d::Tweedie, x::Real; ε::Float64=1e-18, maxiterations::Int=4096)
    #~ Compute estimates of kmax and Wmax
	  kmax = x^(2 - d.p) / (d.ϕ * (d.p - 2))
    α = -shape(d)
    Vtol = log(ε)
    z = (d.p - 1)^α * d.ϕ^(α-1) / (x^α * (d.p - 2))

    #! note: Somehow the approximation does not work at all
    # logVenv(k) = k*(log(z) + (1 - α) - log(k) + α*log(α*k)) + 0.5*log(α)
    logVenv(k) = log(z^k * SpecialFunctions.gamma(1 + α*k) / SpecialFunctions.gamma(1 + k))
    logVmax = logVenv(kmax)
    
    #~ As the terms Vₖ contain a (-1)ᵏ and sinusoidal terms, we work here with the "envelope"
    #  At the end we will re-convert to take the alternating signs into account
    logVs = []
    ks = []
    k = floor(Int, kmax)
    push!(ks, k)
    push!(logVs, logVenv(k))
    
    #~ to the left
    iterations = 1
    while true
        k = ks[end] - 1
        if k < 1 break end
        if (logVmax + Vtol > logVs[end]) || (iterations > maxiterations) break end
        push!(ks, k)
        push!(logVs, logVenv(k))        
    end
    #~ to the right
    k = ceil(Int, kmax)
    push!(ks, k)
    push!(logVs, logVenv(k))
    while true
        k = ks[end] + 1
        if (logVmax + Vtol > logVs[end]) || (iterations > maxiterations) break end
        push!(ks, k)
        push!(logVs, logVenv(k))
    end    

    #~ Compute the Vₖ
    Venv = exp.(logVs)
    signs = [(-1)^(ks[i])*sin(-ks[i]*α*π) for i in eachindex(ks)]
    Vk = Venv .* signs

    #~ Compute the value of the pdf
    a = sum(Vk) / x / π
    P = a * exp((x*θ(d) - κ(d))/d.ϕ)
    return P
end

### Sampling
rand(rng::AbstractRNG, d::Tweedie) = nothing

### Fit model

"""
    fit_mle(::Type{<:Tweedie}, x::AbstractArray{<:Real};
            ϕ0::Real = 1, maxiter::Int = 1000, tol::Real = 1e-16)

Compute the maximum likelihood estimate of the [`Tweedie`](@ref) distribution given fixed
exponent parameter `p` and mean `μ`. Note that estimation of `p` is difficult and a grid search
with error minimization is easier and more stable [see, Dunn & Smyth (2005)]. Additionally, the
mean `μ` can be easily estimated with the sample mean.
"""
function _fit_mle(::Type{<:Tweedie}, x::AbstractArray{<:Real};
    ϕ0::Real = 1, maxiter::Int = 1000, tol::Real = 1e-16)

    N = 0
    nothing
end
