#= Simple module for maximum likelihood estimation
   Has specific models for specific functions as, typically, constraints are different
=#
#/ Start module
module MLEFit

#/ Packages
using Optim
using ForwardDiff
using StatsBase

using Meris

#############################
### SPECIALIZED FUNCTIONS ###

"""
    fitgeneralizedpareto(x::Array{T}; μ=0.0)

Fit a generalized Pareto distribution to data x>=μ [i.e., μ=xmin]
Uses method of moments as the initial guess for the parameters σ and ξ
"""
function fitgeneralizedPareto(x::Array{T}; μ=0.0) where T<:Real
    function negloglikelihood(x, θ)
        logσ, ξ = θ
        σ = exp(logσ)   # ensures σ > 0
        if any(1 .+ ξ .* (x .- μ) ./ σ .<= 0)
            # invalid region
            return Inf
        end
        #~ Return negative log-likelihood
        return -sum(Meris.ParetoLike.loggeneralizedPareto.(x, μ, Ref(σ), Ref(ξ)))
    end
    #~ Use method of moments for initial estimate
    Ex, Vx = mean(x), var(x)
    ξinit = (1 - (Ex - μ)^2 / Vx) / 2
    σinit = (Ex - μ)*(1-ξinit)
    θinit = [log(σinit), ξinit]

    f(θ) = negloglikelihood(x, θ)
    optimres = Optim.optimize(f, θinit, LBFGS(); autodiff=:forward)
    if Optim.converged(optimres)
        σhat, ξhat = optimres.minimizer
        return (; μ=μ, σ=exp(σhat), ξ=ξhat)
    end
    println("Optimizer not converged, returning initial guesses [method of moments]")
    return (; μ=μ, σ=σinit, ξ=ξinit)
end

########################
### HELPER FUNCTIONS ###
function fit_mle(f::Function, data, initial_guess; lower=nothing, upper=nothing, method=Optim.LBFGS())
    if isnothing(lower) || isnothing(upper)
        #~ Unconstrained optimization
        _f = θ -> negloglikelihood(f,θ,data)
	      res = Optim.optimize(_f, initial_guess, method=method, autodiff=:forward)
    else
        #~ Box-constrained optimization
        res = optimize(
            θ -> negloglikelihood(f,θ,data),
            lower, upper,
            initial_guess,
            Fminbox(method),
            autodiff=:forward
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

end # module MLEFit
#/ End module
