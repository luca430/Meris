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

function fit(::Type{TemperedPareto}, x::Array{T}, ε) where {T<:Real}
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
    βinit = 1 / (Ex + ε)
    params = [log(αinit), log(βinit)]    

    #~ Optimize
    #!note: ranges wherein to search are rather arbitrary
    #@TODO: is there a way this can be made more robust?
    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, x),
        [-3.0, log(1/maximum(x))],
        [3.0, log(100/minimum(x))],
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

function fit(::Type{TemperedPareto}, x::Array{T}; εs=nothing) where T<:Real
    xs = sort(x)
    εs = isnothing(εs) ? unique(xs) : εs
    
    αhat = 1.0
    βhat = 1.0
    εhat = 1e-24
    D = Inf
    n = 0
    #/ For each possible xmin in xmins;
    #  - compute the max.-likelihood estimate of the power law exponent γ
    #  - compute the Kolmogorov-Smirnov distance
    #  - extract the xmin for which the MLE γ gives the smallest KS distance
    for i in eachindex(εs)
        ε = εs[i]
        n = count(xs .> ε)
        (n < 50) && (break)        # If less than 50 samples >xmin, break
        #~ Filter data
        _idx = searchsortedfirst(xs, ε) + 1
        _x = xs[_idx:end]
        _P = fit(TemperedPareto, _x, ε)
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        Fv = _ecdf(_x, _x, sorted=true).F     # Values of empirical CDF
        Ftv = max.(0.0, 1.0 .- ccdf.(_P, _x))
        Z = sqrt.(Ftv .* (1 .- Ftv))          # Weight
        distances = abs.(Fv .- Ftv) ./ Z      # Weighted KS distance
        Dhat = maximum(distances)
        #~ If smaller than the current best, update
        if Dhat < D
            αhat = _P.α
            βhat = _P.β
            εhat = _P.ε
            D = Dhat
        end
    end
    return TemperedPareto(αhat, βhat, εhat)
end

########################
### HELPER FUNCTIONS ###
# """
# Compute empirical CDF at points t where F[t] = (no. elements ≤ t) / n
# """
# function _ecdf(xs::Array{T}, t::Array{T}; sorted=false) where T<:Real
#     (!sorted) && (xs = sort(xs))
#     n = length(xs)
#     F = similar(t, Float64)
#     k = 1
#     for i in eachindex(t)
#         #~ Move k until xs[k] > edges[k]
#         while k ≤ n && xs[k] ≤ t[i]
#             k += 1
#         end
#         F[i] = (k-1) / n
#     end
#     return (; F=F, t=t)
# end

# """
# Compute empirical CDF at equally distributed points t
# """
# function _ecdf(x::Array{T}, t::Int; sorted=false) where T<:Real
#     (!sorted) && (xs = sort(x))
#     edges = range(xs[begin], xs[end], length=t) |> collect
#     return _ecdf(x, edges, sorted=true)
# end


# """
# Compute empirical CDF at data points themselves
# """
# function _ecdf(x::Array{T}) where T<:Real
#     return _ecdf(x, sort(x), sorted=false)
# end
