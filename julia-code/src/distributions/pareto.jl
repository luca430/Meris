#= Simple module with formulae for pdfs and cdfs of Pareto distributions
   So far, it includes;
   - Pareto [type I-IV]
     > Pareto Type I --- Pareto
     > Pareto Type II -- Pareto II, and Lomax (for μ=0)
     > Pareto Type III - Pareto III
     > Pareto Type IV -- Pareto IV
   - Generalized Pareto
   - Burr

   Work in progress
   [*] - Tempered Pareto
   [*] - Truncated Pareto
=#
#/ Start module
module ParetoDistribution

using Distributions
using Random
using SpecialFunctions
using StatsBase
using RootSolvers
using Optim

###############
### STRUCTS ###
struct ParetoI{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    ε::T

    ParetoI{T}(α::T, ε::T) where {T<:Real} = new{T}(α, ε)
end

function ParetoI(α::T, ε::T; check_args::Bool = true) where {T<:Real}
	  Distributions.@check_args ParetoI (α, α > 0) (ε, ε > 0)
    return ParetoI{T}(α, ε)
end

struct ParetoII{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    ε::T
    θ::T

    ParetoII{T}(α::T, ε::T, θ::T) where {T<:Real} = new{T}(α, ε, θ)
end

struct ParetoIII{T<:Real} <: ContinuousUnivariateDistribution
    γ::T
    ε::T
    θ::T
    
    ParetoIII{T}(γ::T, ε::T, θ::T) where {T<:Real} = new{T}(γ, ε, θ)
end

struct ParetoIV{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    γ::T
    ε::T
    θ::T

    ParetoIV{T}(α, γ, ε, θ) where {T<:Real} = new{T}(α, γ, ε, θ)
end

function ParetoIV(α::T, γ::T, ε::T, θ::T; check_args::Bool = true) where {T<:Real}
	  Distributions.@check_args ParetoIV (α, α>zero(α)) (γ, γ>zero(γ)) (θ, θ>zero(θ)) (ε,ε>zero(ε))
    return ParetoIV{T}(α,γ,ε,θ)
end

struct Burr{T<:Real} <: ContinuousUnivariateDistribution
    c::T
    α::T
    ε::T

    Burr{T}(c,α,ε) where {T<:Real} = new{T}(c,α,ε)
end

function Burr(c::T, α::T, ε::T; check_args::Bool = true) where {T<:Real}
	  Distributions.@check_args Burr (c, c > zero(c)) (α, α > zero(α)) (ε, ε > zero(ε))
    return Burr{T}(c, α, ε)
end

struct BoundedPareto{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    ε::T
    εmax::T

    BoundedPareto{T}(α,ε,εmax) where {T<:Real} = new{T}(α,ε,εmax)
end

function BoundedPareto(α::T, ε::T, εmax::T; check_args::Bool = true) where {T<:Real}
    Distributions.@check_args BoundedPareto (α, α>zero(α)) (ε,ε>zero(ε)) (εmax,εmax>=ε)
    return BoundedPareto{T}(α,ε,εmax)
end

struct TemperedPareto{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    ε::T

    TemperedPareto{T}(α,β,ε) where {T<:Real} = new{T}(α,β,ε)
end

function TemperedPareto(α::T, β::T, ε::T; check_args::Bool=true) where {T<:Real}
	  Distributions.@check_args TemperedPareto (α, α>zero(α)) (β, β>zero(β)) (ε, ε>zero(ε))
    return TemperedPareto{T}(α, β, ε)
end

"""
    DoublePareto

Continuous double Pareto distribution with switch around typical scale τ, with density function
obtained from the survival function. The density reads

```math
f(x) = [FIX THIS]
```

or, in a more readible format, the density of some random variate `X ~ DoublePareto(α,β,θ,ε)`
[ADD THIS]
"""
struct DoublePareto{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    τ::T
    ε::T

    DoublePareto{T}(α,β,τ,ε) where {T<:Real} = new{T}(α,β,τ,ε)
end

function DoublePareto(α::T, β::T, τ::T, ε::T; check_args::Bool=true) where {T<:Real}
	  Distributions.@check_args DoublePareto (ε, ε>zero(ε)) (α, α>zero(α)) (β, β>=α) (τ,τ>zero(τ))
    return DoublePareto{T}(α,β,τ,ε)
end

"""
    GeneralizedPareto

Restricted form of the generalized Pareto distribution, where the restriction limits it to the
positive domain x≥ε, which means that shape parameter α≥0. In addition, parameters are changed
w.r.t. other sources [like the Wikipedia entry], to remain consistent with the other Pareto
distributions implemented here.

Thus, the density function of the generalized Pareto implemented here reads

```math
f(x) = \\frac{1}{\\theta} \\left(
       1 + \\frac{\\alpha(x - \\varepsilon)}{\\theta}
       \\right)^{-(1+1/\\alpha)}
```

or, in a more readible format, the density of some random variate `X ~ GPD(α,θ,ε)` reads
f(x) = (1/θ)*(1 + α*(x-ε)/θ)^(-(1 + 1/α))

External links
* [Generalized Pareto distribution on Wikipedia](https://en.wikipedia.org/wiki/Generalized_Pareto_distribution)
"""
struct GeneralizedPareto{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    θ::T
    ε::T

    GeneralizedPareto{T}(α,θ,ε) where {T<:Real} = new{T}(α,θ,ε)
end

function GeneralizedPareto(α::T, θ::T, ε::T; check_args::Bool=true) where {T<:Real}
	  Distributions.@check_args GeneralizedPareto (α, α >= zero(α)) (θ, θ >= zero(θ)) (ε, ε >= zero(ε))
    return GeneralizedPareto{T}(α,θ,ε)
end

####################
### CONSTRUCTORS ###

ParetoI(α::Real, ε::Real; check_args::Bool=true) = ParetoI(promote(α,ε)...; check_args=check_args)
ParetoII(α::Real, ε::Real, θ::Real; check_args::Bool=true) = ParetoII(promote(α,ε,θ)...; check_args=check_args)
ParetoIII(γ::Real, ε::Real, θ::Real; check_args::Bool=true) = ParetoIII(promote(γ,ε,θ)...; check_args=check_args)
ParetoIV(α::Real, γ::Real, ε::Real, θ::Real; check_args::Bool=true) = ParetoIV(promote(α,γ,ε,θ)...; check_args=check_args)
TemperedPareto(α::Real, β::Real, ε::Real; check_args::Bool=true) = TemperedPareto(promote(α,β,ε)...; check_args=check_args)
BoundedPareto(α::Real, ε::Real, εmax::Real; check_args::Bool=true) = BoundedPareto(promote(α,ε,εmax)...; check_args=check_args)
GeneralizedPareto(; check_args::Bool=true) = GeneralizedPareto(promote(1.,1.,1.)...; check_args=check_args)
GeneralizedPareto(α::Real, θ::Real, ε::Real; check_args::Bool=true) = GeneralizedPareto(promote(α,θ,ε)...; check_args=check_args)
DoublePareto(α::Real, τ::Real; check_args::Bool=true) = DoublePareto(promote(α,1.,τ,1.)...; check_args=check_args)
DoublePareto(α::Real, β::Real, τ::Real, ε::Real; check_args::Bool=true) = DoublePareto(promote(α,β,τ,ε)...; check_args=check_args)

#~ Burr distribution
Burr(c::Real, α::Real, ε::Real; check_args::Bool=true) = Burr(promote(c,α,ε)...; check_args=check_args)


##################
### STATISTICS ###

mean(d::ParetoI) = if d.α <= 1 return Inf else return d.α*d.ε / (d.α - 1) end
scale(d::ParetoI) = d.ε
shape(d::ParetoI) = d.α
params(d::ParetoI) = (d.α, d.ε)

params(d::BoundedPareto) = (d.α, d.ε, d.εmax)
params(d::DoublePareto) = (d.α, d.β, d.τ, d.ε)


#########################
### DENSITY FUNCTIONS ###
function pdf(d::ParetoI, x::Real)
	  if x < d.ε
        return 0.0
    end
    return d.α*d.ε^d.α * x^(-(d.α+1))
end

function pdf(d::ParetoIV, x::Real)
    if x < d.ε
        return 0.0
    end
    z = (x - d.ε) / d.θ
    return (d.α*d.γ/d.θ) * z^(d.γ-1) * (1 + z^d.γ)^(-1-d.α)
end

function pdf(d::Burr, x::Real)
	  if x <= 0.0
        return 0.0
    end
    return (d.c*d.α / d.ε) * (x / d.ε)^(d.c-1) * (1 + (x / d.ε)^d.c)^(-d.α-1)
end

function pdf(d::TemperedPareto, x::Real)
    if x < d.ε
        return 0.0
    end
	  return d.ε^d.α*exp(d.β*d.ε) * x^(-d.α-1)*exp(-d.β*x) * (d.α + d.β*x)
end

function pdf(d::BoundedPareto, x::Real)
    if x < d.ε || x > d.εmax
        return 0.0
    end
    return d.α*d.ε^d.α * x^(-1-d.α) / (1 - (d.ε/d.εmax)^d.α)
end

function pdf(d::GeneralizedPareto, x::Real)
    if x <= d.ε
        return 0.0
    end
    z = (x - d.ε) / d.θ
    return (1 + d.α*z)^(-1 - 1 / d.α) / d.θ
end

function pdf(d::DoublePareto, x::Real)
    if x < d.ε
        return 0.0
    end
	  εv = (d.ε/d.τ)^(1/d.α) + (d.ε/d.τ)^(1/d.β)
    xα = (x/d.τ)^(1/d.α - 1) / (d.α * d.τ)
    xβ = (x/d.τ)^(1/d.β - 1) / (d.β * d.τ)
    numer = εv * (xα + xβ)
    denum = ((x/d.τ)^(1/d.α) + (x/d.τ)^(1/d.β))^2
    # @info "hm" numer denum
    return numer / denum
end

##########################
### SURVIVAL FUNCTIONS ###
function ccdf(d::ParetoI, x::Real)
	  if x <= d.ε
        return 1.0
    end
    return (d.ε / x)^d.α
end

function ccdf(d::ParetoIV, x::Real)
    if x < d.ε
        return 1.0
    end
    return (1 + ((x - d.ε)/d.θ)^d.γ)^(-d.α)
end

function ccdf(d::Burr, x::Real)
	  if x < d.ε
        return 1.0
    end
    return (1 + (x / d.ε)^d.c)^(-d.α)
end

function ccdf(d::TemperedPareto, x::Real)
	  if x < d.ε
        return 1.0
    end
    return d.ε^d.α*exp(d.β*d.ε) * x^(-d.α)*exp(-d.β*x)
end

function ccdf(d::BoundedPareto, x::Real)
	  (x < d.ε) && (return 1.0)
    (x > d.εmax) && (return 0.0)
    return 1 - (1 - d.ε^d.α*x^(-d.α)) / (1 - (d.ε/d.εmax)^d.α)
end

function ccdf(d::GeneralizedPareto, x::Real)
	  if x < d.ε
        return 1.0
    end
    z = (x - d.ε) / d.θ
    return (1 + d.α*z)^(-1 / d.α)
end

function ccdf(d::DoublePareto, x::Real)
	  return ((d.ε/d.τ)^(1/d.α) + (d.ε/d.τ)^(1/d.β)) / ((x/d.τ)^(1/d.α) + (x/d.τ)^(1/d.β))
end

##########################
### SAMPLING FUNCTIONS ###
"Burr distribution"
xval(d::Burr, u::Real) = d.ε .* ((1 .- u).^(-1/d.α) .- 1) .^ (1 / d.c)
rand(rng::AbstractRNG, d::Burr{T}) where {T<:Real} = xval(d, Random.rand(rng,float(T)))
function rand(rng::AbstractRNG, d::Burr{T}, n::Int) where {T<:Real}
    U = Array{Float64}(undef, n)
    rand!(rng, d, U)
    return U
end
function rand!(rng::AbstractRNG, d::Burr, U::AbstractArray{<:Real})
    Random.rand!(rng, U)
    map!(Base.Fix1(xval, d), U, U)
    return U
end

"Tempered Pareto distribution"
function xval(d::TemperedPareto, u::Real)
    z = log(d.ε^d.α * exp(d.β*d.ε) / (1 - u))
    #~ Define (function, deriviative) for Newton's method
    #! note: we do y = log(x) to avoid Newton's method to have x < 0 [see below]
    fdf(x::Real) = (d.α*x + d.β*exp(x) - z, d.α + d.β*exp(x))
    #~ Solve with good initial guess
    #! note: as x>0, we need the initial guess to be "close", as otherwise the method may
    #        overshoot into x<0 territory, so the guess must depend on `u`.
    guess = log(max(d.ε, z / d.β))
    sol = RootSolvers.find_zero(fdf, RootSolvers.NewtonsMethod{Float64}(guess))
    #~ As we solved for y = log(x), exponentiate
    return exp(sol.root)
end
rand(rng::AbstractRNG, d::TemperedPareto{T}) where {T<:Real} = xval(d, Random.rand(rng,float(T)))
function rand(rng::AbstractRNG, d::TemperedPareto{T}, n::Int) where {T<:Real}
	  U = Array{Float64}(undef, n)
    rand!(rng, d, U)
    return U
end
function rand!(rng::AbstractRNG, d::TemperedPareto{T}, U::AbstractArray{T}) where {T<:Real}
	  Random.rand!(rng, U)
    map!(Base.Fix1(xval, d), U, U)
    return U
end

"Bounded Pareto distribution"
function xval(d::BoundedPareto, u::Real)
	  return (-1 * (u*d.εmax^d.α - u*d.ε^d.α - d.εmax^d.α) / ((d.ε * d.εmax)^d.α))^(-1/d.α)
end
rand(rng::AbstractRNG, d::BoundedPareto{T}) where {T<:Real} = xval(d, Random.rand(rng,float(T)))
function rand(rng::AbstractRNG, d::BoundedPareto{T}, n::Int) where {T<:Real}
	  U = Array{Float64}(undef, n)
    rand!(rng, d, U)
    return U
end
function rand!(rng::AbstractRNG, d::BoundedPareto{T}, U::AbstractArray{T}) where {T<:Real}
	  Random.rand!(rng, U)
    map!(Base.Fix1(xval, d), U, U)
    return U
end

"Double Pareto distribution"
function xval(d::DoublePareto, u::Real)
	  Z = (d.ε/d.τ)^(1/d.α) + (d.ε/d.τ)^(1/d.β)
    #~ Define (function, derivative) for Newton's method
    function fdf(x::Real)
        x = exp(x)
        S = (x / d.τ)^(1/d.α) + (x / d.τ)^(1/d.β)
        f = Z/S - u
        df = -x * (Z / S^2) * ((1/d.α)*(x/d.τ)^(1/d.α) + (1/d.β)*(x/d.τ)^(1/d.β))
        return (f, df)
    end

    guess = log(d.τ * (Z / u)^(d.α))
    sol = RootSolvers.find_zero(fdf, RootSolvers.NewtonsMethod{Float64}(guess))
    return exp(sol.root)
end
rand(rng::AbstractRNG, d::DoublePareto{T}) where {T<:Real} = xval(d, Random.rand(rng,float(T)))
function rand(rng::AbstractRNG, d::DoublePareto{T}, n::Int) where {T<:Real}
	  U = Array{Float64}(undef, n)
    rand!(rng, d, U)
    return U
end
function rand!(rng::AbstractRNG, d::DoublePareto{T}, U::AbstractArray{T}) where {T<:Real}
	  Random.rand!(rng, U)
    map!(Base.Fix1(xval, d), U, U)
    return U
end

#######################
### LOG LIKELIHOODS ###

function logpdf(d::ParetoIV, x::T) where {T<:Real}
    z = (x - d.ε)/d.θ
	  return log(d.α) + log(d.γ) - log(d.θ) + (d.γ-1)*log(z) - (1+d.α)*log(1 + z^d.γ)
end

"Log density function of the tempered Pareto distribution"
function logpdf(d::TemperedPareto, x::T) where {T<:Real}
    return d.α*log(d.ε) + d.β*d.ε - d.β*x - (1+d.α)*log(x) + log(d.α + d.β*x)
end

"Log density function of the bounded Pareto distribution"
function logpdf(d::BoundedPareto, x::T) where {T<:Real}
    if (x < d.ε) || (x > d.εmax)
        return -Inf
    end
	  return (-1-d.α)*log(x) + log(d.α) + d.α*log(d.ε) - log1p(-(d.ε / d.εmax)^d.α)
end



"Log density function of the generalized Pareto distribution"
function logpdf(d::GeneralizedPareto, x::T) where {T<:Real}
    z = (x - d.ε) / d.θ
    (1 + d.α*z <= zero(z)) && (return -Inf)
    #~ Compute logarithm using `log1p` for accuracy
    expn = (-(1 + d.α) / d.α) * log1p(z * d.α)
    return expn - log(d.θ)
end

"Log density function of the double Pareto distribution"
function logpdf(d::DoublePareto, x::T) where {T<:Real}
    (x < d.ε) && (return -Inf)
    #~ No nice analytic form, so just return log of the pdf
    return log(pdf(d, x))
end

###############
### FITTING ###

function fit(::Type{ParetoI}, x::Array{T}; ε=nothing) where {T<:Real}
    xs = sort(x)
    (isnothing(ε)) && (ε = 1.0)
    n = count(x .>= ε)
    (n < 30) && (@warn("Very little data in the tail with ε=$(ε) [only $(n) data points]"))
    #~ Filter data
    idx = searchsortedfirst(xs, ε)
    xfit = xs[idx:end]
    #/ Compute maximum-likelihood estimate of the ParetoI exponent α
    S = sum(log.(xfit / ε))
    αhat = n / S
    return ParetoI(αhat, ε)
end

function fit(::Type{ParetoIV}, x::Array{T}; ε=nothing) where {T<:Real}
	  (isnothing(ε)) && (ε = minimum(x))

    function negloglikelihood(x, params)
        logα, logγ, logθ = params
	      α = exp(logα)
        γ = exp(logγ)
        θ = exp(logθ)
        d = ParetoIV(α, γ, ε, θ)
        return -sum(logpdf.(d, x))
    end

    #~ Init. estimates
    αinit = 1.0
    γinit = 1.0
    θinit = 1.0
    params = [log(αinit), log(γinit), log(θinit)]
    
    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, x),
        params,
        LBFGS(),
        autodiff=:forward
    )
    if Optim.converged(optimres)
        αhat, γhat, θhat = optimres.minimizer
        return ParetoIV(exp(αhat), exp(γhat), ε, exp(θhat))
    end
    @warn("Optimizer not converged, returning initial guesses")
    return ParetoIV(αinit, γinit, ε, θinit)
end

function fit(::Type{TemperedPareto}, x::Array{T}; ε=nothing) where {T<:Real}
    (isnothing(ε)) && (ε = 1.0)

    function negloglikelihood(x, params)
        logα, logβ = params
        α = exp(logα)    #~ ensures α>0
        β = exp(logβ)    #~ ensures β>0
        d = TemperedPareto(α,β,ε)
        return -sum(logpdf.(d, x))
    end

    #~ Initial estimates
    Ex = StatsBase.mean(x)
    αinit = 1 / Ex
    βinit = 1 / (Ex - ε)
    params = [log(αinit), log(βinit)]

    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, x),
        [log(1e-3), log(1e-3)],
        [log(150), log(150)],
        params,
        Fminbox(LBFGS());
        autodiff=:forward
    )
    if Optim.converged(optimres)
        αhat, βhat = optimres.minimizer
        return TemperedPareto(exp(αhat), exp(βhat), ε)
    end
    @warn("Optimizer not converged, returning initial guesses [method of moments]")
    return TemperedPareto(αinit, βinit, ε)
end

function fit(::Type{BoundedPareto}, x::Array{T}; ε=nothing, εmax=nothing) where {T<:Real}
    xs = sort(x)
    (isnothing(ε)) && (ε = minimum(x))
    (isnothing(εmax)) && (εmax = maximum(x))
    #~ Filter data
    minidx = searchsortedfirst(xs, ε)
    maxidx = searchsortedlast(xs, εmax)
    xfit = xs[minidx:maxidx]

    function negloglikelihood(x, params)
        logα = params[begin]
        α = exp(logα)    #~ ensures α>0
        d = BoundedPareto(α,ε,εmax)
        return -sum(logpdf.(d, x))
    end

    #~ Initial estimates as pure Pareto MLE
    S = sum(log.(xfit ./ ε))
    αinit = length(xfit) / S
    if αinit > exp(5)
        return BoundedPareto(αinit, ε, εmax)
    end
    params = [log(αinit)]    

    #~ Use Fminbox as for α→0 the likelihood is extremely flat and infinities will happen
    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, xfit),
        [log(1e-3)],  #~ αmin
        [log(150)],   #~ αmax
        params,
        Fminbox(LBFGS());
        autodiff=:forward
    )
    if Optim.converged(optimres)
        αhat = optimres.minimizer[begin]
        return BoundedPareto(exp(αhat), ε, εmax)
    end
    @warn("Optimizer not converged, returning initial guesses [method of moments]")
    return BoundedPareto(αinit, ε, εmax)
end

function fit(::Type{GeneralizedPareto}, x::Array{T}; ε=nothing) where {T<:Real}
    (isnothing(ε)) && (ε = 1.0)
    
	  function negloglikelihood(x, params)
        logα, logθ = params
        α = 1e-6 + exp(logα)    # ensures α > 0
        θ = 1e-6 + exp(logθ)    # ensures θ > 0
        d = GeneralizedPareto(α, θ, ε)
        if any(1 .+ d.α .* (x .- d.ε) ./ d.θ .<= 0)
            # invalid region
            return Inf
        end
        #~ Return negative log-likelihood
        return -sum(logpdf.(d, x))
    end
    #~ Use method of moments for initial estimate
    Ex, Vx = StatsBase.mean(x), StatsBase.var(x)
    αinit = (1 - (Ex - ε)^2 / Vx) / 2
    θinit = (Ex - ε)*(1-αinit)
    params = [log(αinit), log(θinit)]

    optimres = Optim.optimize(Base.Fix1(negloglikelihood, x), params, LBFGS(); autodiff=:forward)
    if Optim.converged(optimres)
        αhat, θhat = optimres.minimizer
        return GeneralizedPareto(exp(αhat), exp(θhat), ε)
    end
    @warn("Optimizer not converged, returning initial guesses [method of moments]")
    return GeneralizedPareto(αinit, θinit, ε)
end

function fit(::Type{DoublePareto}, x::Array{T}; ε=nothing) where {T<:Real}
    #~ If ε not given, just assume minimum of x
    (isnothing(ε)) && (ε = minimum(x))
    xfit = x[x .> ε]

    function negloglikelihood(x, params)
	      logα, logβ, logτ = params
        α = exp(logα) + 1e-9            #~ Ensures α > 0
        β = α * (1 + exp(logβ) + 1e-9)  #~ Ensures β > α
        τ = exp(logτ)
        d = DoublePareto(α, β, τ, ε)
        return -sum(logpdf.(d, x))
    end
    
    αinit = 0.5
    βinit = 2.0
    τinit = StatsBase.median(x)
    params = [log(αinit), log(βinit/αinit - 1), log(τinit)]

    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, xfit), params,
        LBFGS(), autodiff=:forward
    )
    if Optim.converged(optimres)
        αhat, βhat, τhat = optimres.minimizer
        return DoublePareto(exp(αhat), exp(αhat)*(1 + exp(βhat)), exp(τhat), ε)
    end
    @warn("Optimizer not converged, returning initial guesses")
    return DoublePareto(αinit, βinit, τinit, ε)
end


############
### CDFs ###
function generalizedParetocdf(σ::Float64, ξ::Float64; xmin=0.0)
    function F(x::Float64)
        (x < xmin) && (return 0.0)
        z = (x - xmin) / σ
        if iszero(ξ)
            return 1 - exp(-z)
        end
        return 1 - (1 + ξ*z)^(-1/ξ)
    end
    return F
end

############
### PDFS ###
#~ Generalized Pareto distribution
function generalizedParetopdf(x::Array{T}, σ::Float64, ξ::Float64; xmin=1.0) where T<:Real
	  return generalizedParetopdf.(x, Ref(σ), Ref(ξ); xmin=xmin)
end

function generalizedParetopdf(x::Float64, σ::Float64, ξ::Float64; xmin=1.0)
    (x < xmin) && (return 0.0)
    z = (x - xmin) / σ
    return (1 + ξ*z) ^ (-1 - 1 / ξ) / σ
end





end # module ParetoLike
#/ End module
