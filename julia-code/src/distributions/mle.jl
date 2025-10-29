#= Module to find the parameters of a distribution to closely match the data using MLE
   While for common distributions, such as LogNormal, Gamma, etc., the package `Distributions.jl`
   has the `fit_mle(Function, data)` function, this is not implemented for all distributions.
   Here, a simple MLE scheme is implemented using `Optim.jl` to perform MLE for any pdf f(x,θ)
=#
#/ Start module
module MLEstimator

#/ Packages
using Optim
using SpecialFunctions

#################
### FUNCTIONS ###
function fit(f::Function, data, initial_guess; lower=nothing, upper=nothing, method=Optim.BFGS())
    if isnothing(lower) || isnothing(upper)
        #~ Unconstrained optimization
	      res = Optim.optimize(θ -> negloglikelihood(f,θ,data), initial_guess, method=method)
    else
        #~ Box-constrained optimization
        res = optimize(
            θ -> negloglikelihood(f,θ,data),
            lower, upper,
            initial_guess,
            Fminbox(method)
        )
    end
    θstar = Optim.minimizer(res)
    return θstar
end

########################
### HELPER FUNCTIONS ###
"Define the log-likelihood of a pdf f(x,θ)

 note: assumes that f(x,Ref(θ)) returns density at x for parameter(s) θ
"
function loglikelihood(f::Function, θ, data; eps=1e-32)
    #~ Ensure no infinities appear
    L = log.(clamp.(f.(data, Ref(θ)), eps, 1/eps))
    return sum(L)
end

"Define the negative log-likelihood"
function negloglikelihood(f::Function, θ, data)
    return -loglikelihood(f, θ, data)
end

##########################
### SOME DISTRIBUTIONS ###
"Log-normal distribution"
function LogNormal(x::Float64, θ::Array{Float64})
    μ, σ = θ
    return LogNormal(x, μ, σ)
end

function LogNormal(x::Float64, μ::Float64, σ::Float64)
	  return 1/(x*sqrt(2π*σ^2)) * exp(-(log(x) - μ)^2 / 2 / σ^2)
end

"Beta prime distribution"
function BetaPrime(x::Float64, θ::Array{Float64})
	  α, β = θ
    return BetaPrime(x, exp(α), exp(β))
end

function BetaPrime(x::Float64, α::Float64, β::Float64)
	  return x^(α-1) * (1+x)^(-α-β) / SpecialFunctions.beta(α,β)
end

"Burr distribution"
function Burr(x::Float64, θ::Array{Float64})
    c, k, λ = θ
    return Burr(x, exp(c), exp(k), exp(λ))
end

function Burr(x::Float64, c::Float64, k::Float64, λ::Float64)
	  return c*k/λ * (x/λ)^(c-1) * (1 + (x/λ)^c)^(-k-1)
end

"Shifted Pareto distribution (Lomax distribution)"
function ShiftedPareto(x::Float64, θ::Array{Float64})
    α, λ = θ
    return ShiftedPareto(x, exp(α), exp(λ))
end

function ShiftedPareto(x::Float64, α::Float64, λ::Float64)
    return α*λ^α / (x + λ)^(1+α)
end

"Pareto distribution"
function Pareto(x::Float64, θ::Array{Float64})
    α, ε = θ
    return Pareto(x, α, ε)
end

function Pareto(x::Float64, α::Float64, ε::Float64)
    return x < ε ? 0.0 : α*ε^α / x^(1+α)
end

"Pareto type IV"
function Pareto4(x::Float64, θ::Array{Float64})
    σ, γ, α = θ
    #~ Satisfy some constraints
    return Pareto4(x, exp(σ), exp(γ), exp(α))
end

function Pareto4(x::Float64, σ::Float64, γ::Float64, α::Float64)
    return (α/(σ*γ))*(x/σ)^(1/γ-1) * (1 + (x/σ)^(1/γ))^(-α-1)
    # (x < μ) && (return 0.0)
    # return (α/(γ*σ)) * ((x - μ)/σ)^(1/(γ-1)) * (1 + ((x-μ)/σ)^(1/γ))^(-α-1)
end

"Log-logistic distribution"
function LogLogistic(x::Float64, θ::Array{Float64})
    α, β = θ
    return LogLogistic(x::Float64, exp(α), exp(β))
end

function LogLogistic(x::Float64, α::Float64, β::Float64)
	  return (β/α)*(x/α)^(β-1) / (1 + (x/α)^β)^2
end

"Inverse Gaussian distribution"
function InverseGaussian(x::Float64, θ::Array{Float64})
	  μ, λ = θ
    return InverseGaussian(x, exp(μ), exp(λ))
end

function InverseGaussian(x::Float64, μ::Float64, λ::Float64)
	  return sqrt(λ / (2π * x^3)) * exp(-λ*(x - μ)^2 / (2 * x * μ^2))
end

end # module MLEstimator
#/ End module
