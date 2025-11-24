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

struct TemperedPareto{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    ε::T

    TemperedPareto{T}(α,β,ε) where {T<:Real} = new{T}(α,β,ε)
end

function TemperedPareto(α::T, β::T, ε::T; check_args::Bool=true) where {T<:Real}
	  Distributions.@check_args TemperedPareto (α, α > zero(α)) (β, β > zero(β)) (ε, ε > zero(ε))
    return TemperedPareto{T}(α, β, ε)
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
#~  Pareto distributions
ParetoI(α::Real, ε::Real; check_args::Bool=true) = ParetoI(promote(α,ε)...; check_args=check_args)
ParetoII(α::Real, ε::Real, θ::Real; check_args::Bool=true) = ParetoII(promote(α,ε,θ)..., check_args=check_args)
ParetoIII(γ::Real, ε::Real, θ::Real; check_args::Bool=true) = ParetoIII(promote(γ,ε,θ)..., check_args=check_args)
ParetoIV(α::Real, γ::Real, ε::Real, θ::Real; check_args::Bool=true) = ParetoIV(promote(α,γ,ε,θ)..., check_args=check_args)
TemperedPareto(α::Real, β::Real, ε::Real; check_args::Bool=true) = TemperedPareto(promote(α,β,ε)..., check_args=check_args)
GeneralizedPareto(; check_args::Bool=true) = GeneralizedPareto(promote(1.,1.,1.)...; check_args=check_args)
GeneralizedPareto(α::Real, θ::Real, ε::Real; check_args::Bool=true) = GeneralizedPareto(promote(α,θ,ε)...; check_args=check_args)

#~ Burr distribution
Burr(c::Real, α::Real, ε::Real; check_args::Bool=true) = Burr(promote(c,α,ε)...; check_args=check_args)


##################
### STATISTICS ###

mean(d::ParetoI) = if d.α <= 1 return Inf else return d.α*d.ε / (d.α - 1) end
scale(d::ParetoI) = d.ε
shape(d::ParetoI) = d.α
params(d::ParetoI) = (d.α, d.ε)

#########################
### DENSITY FUNCTIONS ###
function pdf(d::ParetoI, x::Real)
	  if x <= d.ε
        return 0.0
    end
    return d.α*d.ε^d.α * x^(-(d.α+1))
end

function pdf(d::Burr, x::Real)
	  if x <= d.ε
        return 0.0
    end
    return (d.c*d.α / d.ε) * (x / d.ε)^(d.c-1) * (1 + (x / d.ε)^d.c)^(-d.α-1)
end

function pdf(d::TemperedPareto, x::Real)
    if x <= d.ε
        return 0.0
    end
	  return d.ε^d.α*exp(d.β*d.ε) * x^(-d.α-1)*exp(-d.β*x) * (d.α + d.β*x)
end

function pdf(d::GeneralizedPareto, x::Real)
    if x <= d.ε
        return 0.0
    end
    z = (x - d.ε) / d.θ
    return (1 + d.α*z)^(-1 - 1 / d.α) / d.θ
end

##########################
### SURVIVAL FUNCTIONS ###
function ccdf(d::ParetoI, x::Real)
	  if x <= d.ε
        return 1.0
    end
    return (d.ε / x)^d.α
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

function ccdf(d::GeneralizedPareto, x::Real)
	  if x < d.ε
        return 1.0
    end
    z = (x - d.ε) / d.θ
    return (1 + d.α*z)^(-1 / d.α)
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
    #~ Define (function, deriviative) for NewtonsMethod
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

#######################
### LOG LIKELIHOODS ###

"Log density function of the tempered Pareto distribution"
function logpdf(d::TemperedPareto, x::T) where {T<:Real}
    return d.α*log(d.ε) + d.β*d.ε - d.β*x - (1+d.α)*log(x) + log(d.α + d.β*x)
end


"Log density function of the generalized Pareto distribution, with an expansion near zero."
function logpdf(d::GeneralizedPareto, x::T) where T<:Real
    z = (x - d.ε) / d.θ
    (1 + d.α*z <= zero(z)) && (return -Inf)
    #~ Compute logarithm using `log1p` for accuracy
    expn = (-(1 + d.α) / d.α) * log1p(z * d.α)
    return expn - log(d.θ)
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
    αhat = 1 + n / S
    return ParetoI(αhat, ε)
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

    optimres = Optim.optimize(Base.Fix1(negloglikelihood, x), params, LBFGS(); autodiff=:forward)
    if Optim.converged(optimres)
        αhat, βhat = optimres.minimizer
        return TemperedPareto(αhat, βhat, ε)
    end
    @warn("Optimizer not converged, returning initial guesses [method of moments]")
    return TemperedPareto(αinit, βinit, ε)
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
