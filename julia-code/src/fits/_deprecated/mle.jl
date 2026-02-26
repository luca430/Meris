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
    fittemperedPareto(x::Array{T}; xmin=1.0)

Fit a tempered Pareto distribution to data x≥xmin
"""
function fittemperedPareto(x::Array{T}; xmin=1.0) where {T <: Real}
	  function negloglikelihood(x, θ)
        logα, logβ = θ
        α, β = exp(logα), exp(logβ)    #~ ensures α,β > 0
        return -sum(Meris.ParetoLike.logtemperedPareto.(x, Ref(α), Ref(β), Ref(xmin)))
    end
    #~ Initial estimates
    Ex = mean(x)
    αinit = 1 / Ex
    βinit = 1 / (Ex - xmin)
    θinit = [log(αinit), log(βinit)]

    f(θ) = negloglikelihood(x, θ)
    optimres = Optim.optimize(f, θinit, LBFGS(); autodiff=:forward)
    if Optim.converged(optimres)        
        αhat, βhat = optimres.minimizer
        return (; α=exp(αhat), β=exp(βhat))
    end
    println("Optimizer not converged, return initial guesses")
    return (; α=αinit, β=βinit)
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
