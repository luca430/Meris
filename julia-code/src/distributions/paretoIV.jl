###############
### STRUCTS ###

struct ParetoIV{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    θ::T
    ε::T

    ParetoIV{T}(α,β,θ,ε) where {T<:Real} = new{T}(α,β,θ,ε)
end

function ParetoIV(α::T, β::T, θ::T, ε::T; check_args::Bool = true) where {T<:Real}
	  Distributions.@check_args ParetoIV (α, α>zero(α)) (β, β>zero(β)) (θ, θ>zero(θ)) (ε,ε>=zero(ε))
    return ParetoIV{T}(α,β,θ,ε)
end

####################
### CONSTRUCTORS ###

ParetoII(α::Real, θ::Real, ε::Real; check_args::Bool=true) = ParetoIV(promote(α,1.0,θ,ε)...; check_args=check_args)
Lomax(α::Real, θ::Real; check_args::Bool=true) = ParetoIV(promote(α,1.0,θ,0.0)...; check_args=check_args)
ParetoIII(β::Real, θ::Real, ε::Real; check_args::Bool=true) = ParetoIV(promote(1.0,β,θ,ε)...; check_args=check_args)
Burr(α::Real, β::Real, θ::Real; check_args::Bool=true) = ParetoIV(promote(α,β,θ,0.0)...; check_args=check_args)
ParetoIV(α::Real, β::Real, θ::Real, ε::Real; check_args::Bool=true) = ParetoIV(promote(α,β,θ,ε)...; check_args=check_args)

##################
### STATISTICS ###

params(d::ParetoIV) = (d.α, d.β, d.θ, d.ε)

#########################
### DENSITY FUNCTIONS ###

function pdf(d::ParetoIV, x::Real)
    if x < d.ε
        return 0.0
    end
    z = (x - d.ε) / d.θ
    return (d.α / (d.β * d.θ)) * z^(1/d.β - 1) * (1 + z^(1/d.β))^(-1-d.α)
end

function logpdf(d::ParetoIV, x::T) where {T<:Real}
    (x < d.ε) && (return -Inf)
    z = (x - d.ε) / d.θ
	  return log(d.α) - log(d.β) - log(d.θ) + (1/d.β - 1)*log(z) - (1 + d.α)*log1p(z^(1/d.β))
end

#########################
### SURVIVAL FUNCTION ###

function ccdf(d::ParetoIV, x::Real)
	  (x < d.ε) && (return 1.0)
    z = (x - d.ε)/d.θ
    return (1 + z^(1/d.β))^(-d.α)
end

##########################
### SAMPLING FUNCTIONS ###
xval(d::ParetoIV, u::Real) = d.ε + d.θ * ((1 - u)^(-1/d.α) - 1)^(d.β)
rand(rng::AbstractRNG, d::ParetoIV{T}) where {T<:Real} = xval(d, Random.rand(rng,float(T)))
function rand(rng::AbstractRNG, d::ParetoIV{T}, n::Int) where {T<:Real}
    return map(Base.Fix1(xval, d), Random.rand(rng,float(T),n))    
end
function rand!(rng::AbstractRNG, d::ParetoIV, U::AbstractArray{<:Real})
    Random.rand!(rng, U)
    map!(Base.Fix1(xval, d), U, U)
    return U
end

###############
### FITTING ###

"""
    @TODO Ensure some parameters remain fixed when subtypes are called
          Most likely, we need seperate structures for this anyways, and then
          write `fit` functions for each.
"""
function fit(::Type{ParetoIV}, x::Array{T}; ε=nothing) where {T<:Real}
	  (isnothing(ε)) && (ε = minimum(x))

    function negloglikelihood(x, params)
        logα, logβ, logθ = params
	      α = exp(logα)
        β = exp(logβ)
        θ = exp(logθ)
        d = ParetoIV(α, β, θ, ε)
        return -sum(logpdf.(d, x))
    end

    #/ Init. estimates
    #~ Quantile estimator for θ
    #  θ determines the "distance" from ε to the start of the heavy tail.
    #  So taking a 10% quantile seems a decent educated guess.
    θinit = quantile(x .- ε, 0.1)
    #~ Hill estimator for α for the top 16% of values
    αinit = 1.0 / hills_estimator(x, sorted=false)[end]
    #~ Guess of β, typically β∈[0.5,1.0]
    βinit = 0.75

    params = [log(αinit), log(βinit), log(θinit)]

    #~ Optimize
    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, x),
        [log(1e-3), log(1e-3), log(minimum(x) .- ε)],
        [log(10.0), log(10.0), log(maximum(x))],
        params,
        Fminbox(LBFGS()),
        autodiff=:forward
    )
    if Optim.converged(optimres)
        αhat, βhat, θhat = optimres.minimizer
        return ParetoIV(exp(αhat), exp(βhat), exp(θhat), ε)
    end
    @warn("Optimizer not converged, returning initial guesses")
    return ParetoIV(αinit, βinit, θinit, ε)
end

########################
### HELPER FUNCTIONS ###

"""
    hills_estimator

Naive Hill's estimator for data
When `sorted=true`, expects the data to be sorted [from large to small!]
"""
function hills_estimator(x::Array{T}; sorted=false) where T<:Real
	  xs = sorted ? x : reverse(sort(x))
    S = log.(xs)
    C = cumsum(S)
    ξ = similar(xs, Float64)
    for k in eachindex(xs)
        (k == 1) && (continue)
        ξ[k] = (C[k-1] - (k-1)*S[k]) / (k - 1)
    end
    return ξ
end
