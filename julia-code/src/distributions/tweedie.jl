###############
### STRUCTS ###

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
    b::T
    μ::T
    ϕ::T

    Tweedie{T}(b::T, μ::T, ϕ::T) where {T<:Real} = new{T}(b, μ, ϕ)
end

function Tweedie(b::T, μ::T, ϕ::T; check_args::Bool=true) where {T<:Real}
    #~ Special cases
    (b == 0) && (return Normal(μ, sqrt(ϕ)))
    (b == 1) && (return Poisson(μ))
    (b == 2) && (return Gamma(1/ϕ,ϕ*μ))
    # (b == 3) && (return InverseGaussian(μ,1/ϕ))
	  Distributions.@check_args Tweedie (b, !(zero(b)<b<one(b))) (μ, μ>=zero(μ)) (ϕ, ϕ>zero(ϕ))
    return Tweedie{T}(b, μ, ϕ)
end

### Constructors
Tweedie(b::Real; check_args::Bool=true) = Tweedie(b, one(b), one(b); check_args=check_args)
Tweedie(b::Real, μ::Real, ϕ::Real; check_args::Bool=true) = Tweedie(promote(b, μ, ϕ)...; check_args=check_args)

### Parameters
params(d::Tweedie) = (d.b, d.μ, d.ϕ)

### Statistics
mean(d::Tweedie) = d.μ
var(d::Tweedie) = d.μ^d.b
shape(d::Tweedie) = (2 - d.b) / (d.b - 1)
θ(d::Tweedie) = isone(d.b) ? log(d.μ) : d.μ^(1 - d.b) / (1 - d.b)
κ(d::Tweedie) = isone(d.b - 1) ? log(d.μ) : d.μ^(2 - d.b) / (2 - d.b)

#########################
### DENSITY FUNCTIONS ###

"""
    pdf(d::Tweedie, x::Real)

The density function of a Tweedie distribution generally does not have a closed form (depending
on the value of `p`), and hence the density at `x` needs to be approximated. Here, the series
expansion of Dunn & Smyth (2005) is used. In general, the series expansion is good for values
`x` not too small (a rough estimate is x>0.1) and is excellent for x≥0.5. Simultaneously, the
series expansion, and also the maximum likelihood estimate, becomes worse as `p` approaches the
Poisson limit at p=1.
"""

function pdf(d::Tweedie, x::Real; ε::Float64 = 1e-16, maxiterations::Int=8192)
	  if iszero(x)
        return exp(-(d.μ^(2 - d.b)) / (d.ϕ * (2 - d.b)))
    end
    if d.b > 2
        return _pdfpoissonstable(d::Tweedie, x::Real; ε=ε, maxiterations=maxiterations)
    end
    return _pdfpoissongamma(d::Tweedie, x::Real; ε=ε, maxiterations=maxiterations)
end

"""
Approximate the density for a Tweedie distribution evaluated at `x`
"""
function _pdfpoissongamma(d::Tweedie, x::Real; ε::Float64=1e-16, maxiterations::Int=4096)
    _, logWs = _seriesexpansion1b2(d, x; ε=ε, maxiterations=maxiterations)
    #~ Compute the Wₖ
    Wk = exp.(logWs)
    #~ Compute the value of the pdf
    a = sum(Wk) / x
    return a * exp((x*θ(d) - κ(d)) / d.ϕ)
end

function _pdfpoissonstable(d::Tweedie, x::Real; ε::Float64=1e-16, maxiterations::Int=4096)
    ks, Vs = _seriesexpansionb2(d, x; ε=ε, maxiterations=maxiterations)
    #~ Compute the Vₖ
    angles = mod.(ks*shape(d)*π, 2π)
    signs = [(-1)^ks[i] * sin(angles[i]) for i in eachindex(ks)]
    Vk = Vs .* signs
    #~ Compute the value of the pdf
    a = sum(Vk) / x / π
    return a * exp((x*θ(d) - κ(d))/d.ϕ)
end

"""
Approximate the infinite series W(x;ϕ,p) = ∑ₖ Wₖ using methods from Dunn & Smyth (2005)
Specialized version for 1<b<2
"""
function _seriesexpansion1b2(d::Tweedie, x::Real; ε::Float64=1e-16, maxiterations::Int=4096)
    #~ Compute estimates of kmax and Wmax
    kmax = x^(2 - d.b) / (d.ϕ * (2 - d.b))
    α = -shape(d)
    Wtol = log(ε)
    z = x^(-α) * (d.b - 1)^α / (d.ϕ^(1 - α) * (2 - d.b))

    logWk(k) = k > 32 ?
               log(z^k / SpecialFunctions.gamma(-k*α) / SpecialFunctions.gamma(1+k)) :
               k*(log(z) + (1-α) + α*log(-α) - (1-α)*log(k)) - log(2π) - 0.5*log(-α) - log(k)
    logWmax = logWk(kmax)

    #~ Instantiate and allocate
    iterations = 0
    logWs = []
    ks = []
    k = floor(Int, kmax)
    
    if k > 0
        push!(ks, k)
        push!(logWs, logWk(k))
    
        #~ to the left        
        while true
            k = ks[end] - 1
            if k < 1 break end
            iterations += 1
            if (logWmax + Wtol > logWs[end]) || (iterations > maxiterations) break end
            push!(ks, k)
            push!(logWs, logWk(k))
        end
    end

    #~ to the right
    k = ceil(Int, kmax)
    push!(ks, k)
    push!(logWs, logWk(k))
    while true
        k = ks[end] + 1
        iterations += 1
        if (logWmax + Wtol > logWs[end]) || (iterations > maxiterations) break end
        push!(ks, k)
        push!(logWs, logWk(k))
    end

    (iterations > maxiterations) && (@warn("Sum did not converge, consider raising iterations."))

    return (ks, logWs)
end

"""
Approximate the infinite series W(x;ϕ,p) = ∑ₖ Wₖ using methods from Dunn & Smyth (2005)
Specialized version for b>2
"""
function _seriesexpansionb2(d::Tweedie, x::Real; ε::Float64=1e-18, maxiterations::Int=4096)
    #~ Compute estimates of kmax and Wmax
	  kmax = x^(2 - d.b) / (d.ϕ * (d.b - 2))
    α = -shape(d)
    Vtol = log(ε)
    z = (d.b - 1)^α * d.ϕ^(α-1) / (x^α * (d.b - 2))

    #~ Define the envelope function
    #! note: for large k, use the approximation instead, yet for small k use the true value
    logVenv(k) = k > 32 ?
                 k*(log(z) + (1 - α) - log(k) + α*log(α*k)) + 0.5*log(α) :
                 log(z^k * SpecialFunctions.gamma(1 + α*k) / SpecialFunctions.gamma(1 + k))
    logVmax = logVenv(kmax)
    
    #~ As the terms Vₖ contain a (-1)ᵏ and sinusoidal terms, we work here with the "envelope"
    #  At the end we will re-convert to take the alternating signs into account
    logVs = []
    ks = []
    k = floor(Int, kmax)
    push!(ks, k)
    push!(logVs, logVenv(k) - logVmax)
    
    #~ to the left
    iterations = 0
    while true
        k = ks[end] - 1
        if k < 1 break end
        iterations += 1
        if (logVs[end] < Vtol) || (iterations > maxiterations) break end
        push!(ks, k)
        push!(logVs, logVenv(k) - logVmax)
    end
    #~ to the right
    k = ceil(Int, kmax)
    push!(ks, k)
    push!(logVs, logVenv(k) - logVmax)
    while true
        k = ks[end] + 1
        iterations += 1
        if (logVs[end] < Vtol) || (iterations > maxiterations) break end
        push!(ks, k)
        push!(logVs, logVenv(k) - logVmax)
    end
    
    (iterations > maxiterations) && (@warn("Sum did not converge, consider raising iterations."))
    #~ Compute logV
    Vs = exp.(logVs) * exp(logVmax)
    return (ks, Vs)
end

function logpdf(d::Tweedie, x; ε::Float64=1e-18, maxiterations::Int=4096)
    return log(pdf(d, x; ε=ε, maxiterations=maxiterations))
end

### Sampling
rand(rng::AbstractRNG, d::Tweedie) = nothing


### Fit model

"""
    fit_mle(::Tweedie, x::AbstractArray{<:Real};
            ϕ0::Real = 1, maxiter::Int = 1000, tol::Real = 1e-16)

Compute the maximum likelihood estimate of the [`Tweedie`](@ref) distribution given fixed
exponent parameter `b` and mean `μ`. Note that estimation of `b` is difficult and a grid search
with error minimization is easier and more stable [see, Dunn & Smyth (2005)]. Additionally, the
exponent `b` can be obtained more robustly from Taylor's law itself. That is, if Taylor's law is
not been established, then fitting a Tweedie distribution is probably not a good idea. However,
if it is observed and the exponent `b` has been obtained, the MLE for ϕ can be obtained.
Finally note that the mean `μ` can of course be easily estimated with the sample mean.
"""
function fit(
    d::Tweedie, x::AbstractArray{<:Real};
    ϕ0::Float64 = StatsBase.var(x) / (StatsBase.mean(x)^d.b),
    maxiterations::Int = 4096,
    ε::Real = 1e-16
)
    _f(d,x) = logpdf(d, x; ε=ε, maxiterations=maxiterations)
    _df(d,x) = _difflogpdf(d, x; ε=ε, maxiterations=maxiterations)
    uf = similar(x)
    udf = similar(x)
    function f(logϕ)
        Tw = Tweedie(d.b, d.μ, exp(logϕ))
        return sum(map!(Base.Fix1(_f, Tw), uf, x)) * exp(logϕ)
    end
    function df(logϕ)
        Tw = Tweedie(d.b, d.μ, exp(logϕ))
        return sum(map!(Base.Fix1(_df, Tw), udf, x)) * exp(logϕ)
    end
    sol = RootSolvers.find_zero(logϕ -> (f(logϕ), df(logϕ)), NewtonsMethod{Float64}(log(ϕ0)))
    return sol
end

"[private] Derivative of the log density of a Tweedie distribution"
function _difflogpdf(d::Tweedie, x::Real; ε::Float64=1e-18, maxiterations::Int=4096)
    if iszero(x)
        return d.μ^(2-d.b) / (d.ϕ^2 * (2 - d.b))
    end
    if d.b > 2
        return _difflogb2(d, x; ε=ε, maxiterations=maxiterations)
    end
    return nothing
    # return _diffseriesexpansion1b2(d, x; ε=ε, maxiterations=maxiterations)
end

function _difflogb2(d::Tweedie, x::Real; ε::Float64=1e-18, maxiterations::Int=4096)
	  ks, Vs = _seriesexpansionb2(d, x; ε=ε, maxiterations=maxiterations)
    #~ Compute the Vₖ
    angles = mod.(ks*shape(d)*π, 2π)
    signs = [(-1)^ks[i] * sin(angles[i]) for i in eachindex(ks)]
    Vk = Vs .* signs
    #~ Compute the value of the derivative of the log pdf
    V = (-1 - shape(d)) / d.ϕ
    V = V * sum(ks .* Vk) / sum(Vk)
    return (x * d.μ^(1-d.b)) / (d.ϕ^2 * (d.b - 1)) + d.μ^(2 - d.b) / (d.ϕ^2 * (2 - d.b)) + V
end

"Series expansion for the derivative of the density w.r.t. ϕ"
function _diffseriesexpansionb2(d::Tweedie, x::Real; ε::Float64=1e-18, maxiterations::Int=4096)
    #~ Compute estimates of kmax and Wmax
	  kmax = x^(2 - d.b) / (d.ϕ * (d.b - 2))
    α = -shape(d)
    Vtol = log(ε)
    z = (d.b - 1)^α * d.ϕ^(α-1) / (x^α * (d.b - 2))

    #~ Define the envelope function
    #! note: for large k, use the approximation instead, yet for small k use the true value
    logVenv(k) = k > 24 ?
                 k*(log(z) + (1 - α) - log(k) + α*log(α*k)) + 0.5*log(α) :
                 log(z^k * SpecialFunctions.gamma(1 + α*k) / SpecialFunctions.gamma(1 + k))
    logkVmax = logVenv(kmax) + log(kmax)
    
    #~ As the terms Vₖ contain a (-1)ᵏ and sinusoidal terms, we work here with the "envelope"
    #  At the end we will re-convert to take the alternating signs into account
    logVs = []
    ks = []
    k = floor(Int, kmax)
    k = k < 1 ? 1 : k    #~ check if kmax<1, then start at 1
    push!(ks, k)
    push!(logVs, logVenv(k) + log(k) - logkVmax)
    #~ to the left
    iterations = 0
    while k > 1
        k = ks[end] - 1
        # if k < 1 break end
        iterations += 1
        if (logVs[end] < Vtol) || (iterations > maxiterations) break end
        push!(ks, k)
        push!(logVs, logVenv(k) + log(k) - logkVmax)
    end
    #~ to the right
    k = ceil(Int, kmax)
    k = (k in ks) ? k + 1 : k    #~ if kmax < 1, start at kmax+1
    push!(ks, k)
    push!(logVs, logVenv(k) + log(k) - logkVmax)
    while true
        k = ks[end] + 1
        iterations += 1
        if (logVs[end] < Vtol) || (iterations > maxiterations) break end
        push!(ks, k)
        push!(logVs, logVenv(k) + log(k) - logkVmax)
    end
    
    (iterations > maxiterations) && (@warn("Sum did not converge, consider raising iterations."))
    #~ Compute logV
    Vs = exp.(logVs) * exp(logkVmax)
    return (ks, Vs)
end
