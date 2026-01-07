struct TemperedPareto{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    ε::T

    TemperedPareto{T}(α,β,ε) where {T<:Real} = new{T}(α,β,ε)
end

function TemperedPareto(α::T, β::T, ε::T; check_args::Bool=true) where {T<:Real}
	  @check_args TemperedPareto (α, α>zero(α)) (β, β>zero(β)) (ε, ε>zero(ε))
    return TemperedPareto{T}(α, β, ε)
end

####################
### CONSTRUCTORS ###

TemperedPareto(α::Real, β::Real, ε::Real; check_args::Bool=true) = TemperedPareto(promote(α,β,ε)...; check_args=check_args)

##################
### STATISTICS ###
params(d::TemperedPareto) = (d.α, d.β, d.ε)

#########################
### DENSITY FUNCTIONS ###

function pdf(d::TemperedPareto, x::Real)
    if x < d.ε
        return 0.0
    end
	  return d.ε^d.α*exp(d.β*d.ε) * x^(-d.α-1)*exp(-d.β*x) * (d.α + d.β*x)
end

function logpdf(d::TemperedPareto, x::T) where {T<:Real}
    return d.α*log(d.ε) + d.β*d.ε - d.β*x - (1+d.α)*log(x) + log(d.α + d.β*x)
end

##########################
### SURVIVAL FUNCTIONS ###

function ccdf(d::TemperedPareto, x::Real)
	  if x < d.ε
        return 1.0
    end
    return d.ε^d.α*exp(d.β*d.ε) * x^(-d.α)*exp(-d.β*x)
end

##########################
### SAMPLING FUNCTIONS ###

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
    return map(Base.Fix1(xval, d), Random.rand(rng,float(T),n))    
end
function rand!(rng::AbstractRNG, d::TemperedPareto, U::AbstractArray{<:Real})
    Random.rand!(rng, U)
    map!(Base.Fix1(xval, d), U, U)
    return U
end

###############
### FITTING ###

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
    #  uses standard Pareto MLE for α
    Ex = StatsBase.mean(x)
    xs = sort(x)
    idx = searchsortedfirst(xs, ε)
    xfit = xs[idx:end]
    S = sum(log.(xfit / ε))
    αinit = length(x) / S
    βinit = 1 / (Ex - ε)
    params = [log(αinit), log(βinit)]    

    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, x),
        [-3.0, log(1/maximum(x))],
        [3.0, 0.0],
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
