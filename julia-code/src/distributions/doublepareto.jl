###############
### STRUCTS ###

struct DoublePareto{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    τ::T
    ε::T

    DoublePareto{T}(α::T, β::T, τ::T, ε::T) where {T<:Real} = new{T}(α,β,τ,ε)
end

function DoublePareto(α::T, β::T, τ::T, ε::T; check_args::Bool = true) where {T<:Real}
	  @check_args DoublePareto (α, α>zero(α)) (α,α>β) (ε, ε>zero(ε)) (τ, τ>ε)
    return DoublePareto{T}(α,β,τ,ε)
end

####################
### CONSTRUCTORS ###

DoublePareto(α::Real,β::Real,τ::Real,ε::Real; check_args::Bool=true) =
    DoublePareto(promote(α,β,τ,ε)...; check_args=check_args)


##################
### STATISTICS ###

params(d::DoublePareto) = (d.α, d.β, d.τ, d.ε)

########################
### DENSITY FUNCTION ###

function pdf(d::DoublePareto, x::Real)
    (x < d.ε) && (return 0.0)
    T = d.τ^(-1-d.α) / d.τ^(-1-d.β)
    Z = d.τ^(-d.α) / d.α + T * (d.ε^(-d.β) - d.τ^(-d.β)) / d.β
    
    A = T / Z
    B = 1 / Z
    return ifelse(x <= d.τ, A*x^(-1-d.β), B*x^(-1-d.α))
end

############################
### LOG DENSITY FUNCTION ###
function logpdf(d::DoublePareto, x::Real)
	  (x < d.ε) && (return -Inf)
    T = d.τ^(-1-d.α) / d.τ^(-1-d.β)
    Z = d.τ^(-d.α) / d.α + T * (d.ε^(-d.β) - d.τ^(-d.β)) / d.β
    
    A = T / Z
    B = 1 / Z
    return ifelse(x <= d.τ, log.(A) - (1 + d.β)*log(x), log(B) - (1 + d.α)*log(x))
end

function logpdf(d::DoublePareto, x::Vector)
    T = d.τ^(-1-d.α) / d.τ^(-1-d.β)
    Z = d.τ^(-d.α) / d.α + T * (d.ε^(-d.β) - d.τ^(-d.β)) / d.β
    
    logA = log(T / Z)
    logB = log(1 / Z)

    return ifelse.(
        x .< d.ε, -Inf,
        ifelse.(
            x .<= d.τ,
            logA .- (1 + d.β).*log.(x),
            logB .- (1 + d.α).*log.(x)
        )
    )
end

#########################
### SURVIVAL FUNCTION ###

function ccdf(d::DoublePareto, x::Real)
	  (x < d.ε) && (return 1.0)
    T = d.τ^(-1-d.α) / d.τ^(-1-d.β)
    Z = d.τ^(-d.α) / d.α + T * (d.ε^(-d.β) - d.τ^(-d.β)) / d.β
    
    A = T / Z
    B = 1 / Z
    return 1 - ifelse(x <= d.τ,
        (A / d.β) * (d.ε^(-d.β) - x^(-d.β)),
        (A / d.β) * (d.ε^(-d.β) - d.τ^(-d.β)) + (B / d.α)*(d.τ^(-d.α) - x^(-d.α))
    )
end

################
### SAMPLING ###

###############
### FITTING ###

function fit(::Type{DoublePareto}, x::Array{T}; ε=nothing) where {T<:Real}
    (isnothing(ε)) && (ε = minimum(x))

    function negloglikelihood(x, params)
	      logα, logβ, logτ = params
        α = exp(logα)
        β = α - exp(logβ)
        τ = ε + exp(logτ)
        d = DoublePareto(α, β, τ, ε)
        return -sum(logpdf(d,x))
    end

    #~ Initial estimates
    αinit = 1.1
    βinit = 0.5
    τinit = median(x)
    params = [log(αinit), log(βinit), log(τinit)]

    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, x),
        [log(1e-3), log(1e-3), log(minimum(x))],
        [log(5.0), log(5.0), log(maximum(x))],
        params,
        Fminbox(LBFGS()),
        autodiff=:forward
    )
    if Optim.converged(optimres)
        αopt, βopt, τopt = optimres.minimizer
        αhat = exp(αopt)
        βhat = αhat - exp(βopt)
        τhat = ε + exp(τopt)
        return DoublePareto(αhat, βhat, τhat, ε)
    end
    @warn("Optimizer not converged, returning initial guesses")
    return DoublePareto(αinit, βinit, τinit, ε)
end
