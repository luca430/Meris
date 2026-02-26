#= GeneralizedPareto;
   most of the code is yoinked from Distributions.jl, but ours includes a fitting procedure
=#

###############
### STRUCTS ###
"""
    GeneralizedPareto(ε, σ, ξ)

The *Generalized Pareto distribution* (GPD) with shape parameter `ξ`, scale `σ` and location `μ` has probability density function

```math
f(x; \\varepsilon	, \\sigma, \\xi) = \\begin{cases}
        \\frac{1}{\\sigma}(1 + \\xi \\frac{x - \\varepsilon	}{\\sigma} )^{-\\frac{1}{\\xi} - 1} & \\text{for } \\xi \\neq 0 \\\\
        \\frac{1}{\\sigma} e^{-\\frac{\\left( x - \\varepsilon	 \\right) }{\\sigma}} & \\text{for } \\xi = 0
    \\end{cases}~,
    \\quad x \\in \\begin{cases}
        \\left[ \\varepsilon	, \\infty \\right] & \\text{for } \\xi \\geq 0 \\\\
        \\left[ \\varepsilon	, \\varepsilon	 - \\sigma / \\xi \\right] & \\text{for } \\xi < 0
    \\end{cases}
```

```julia
GeneralizedPareto()             # GPD with unit shape and unit scale, i.e. GeneralizedPareto(0, 1, 1)
GeneralizedPareto(ξ)            # GPD with shape ξ and unit scale, i.e. GeneralizedPareto(0, 1, ξ)
GeneralizedPareto(σ, ξ)         # GPD with shape ξ and scale σ, i.e. GeneralizedPareto(0, σ, ξ)
GeneralizedPareto(μ, σ, ξ)      # GPD with shape ξ, scale σ and location μ.

params(d)       # Get the parameters, i.e. (ε, σ, ξ)
location(d)     # Get the location parameter, i.e. ε
scale(d)        # Get the scale parameter, i.e. σ
shape(d)        # Get the shape parameter, i.e. ξ
```

External links

* [Generalized Pareto distribution on Wikipedia](https://en.wikipedia.org/wiki/Generalized_Pareto_distribution)

"""

struct GeneralizedPareto{T<:Real} <: ContinuousUnivariateDistribution
    ε::T
    σ::T
    ξ::T
    GeneralizedPareto{T}(ε::T, σ::T, ξ::T) where {T} = new{T}(ε, σ, ξ)
end

function GeneralizedPareto(ε::T, σ::T, ξ::T; check_args::Bool=true) where {T <: Real}
    @check_args GeneralizedPareto (σ, σ > zero(σ))
    return GeneralizedPareto{T}(ε, σ, ξ)
end

####################
### CONSTRUCTORS ###
function GeneralizedPareto(ε::Real, σ::Real, ξ::Real; check_args::Bool=true)
    return GeneralizedPareto(promote(ε, σ, ξ)...; check_args=check_args)
end

function GeneralizedPareto(ε::Integer, σ::Integer, ξ::Integer; check_args::Bool=true)
    GeneralizedPareto(float(ε), float(σ), float(ξ); check_args=check_args)
end

function GeneralizedPareto(σ::Real, ξ::Real; check_args::Bool=true)
    GeneralizedPareto(zero(σ), σ, ξ; check_args=check_args)
end

function GeneralizedPareto(ξ::Real; check_args::Bool=true)
    GeneralizedPareto(zero(ξ), one(ξ), ξ; check_args=check_args)
end

GeneralizedPareto() = GeneralizedPareto{Float64}(0.0, 1.0, 1.0)

##################
### Parameters ###
location(d::GeneralizedPareto) = d.ε
scale(d::GeneralizedPareto) = d.σ
shape(d::GeneralizedPareto) = d.ξ
params(d::GeneralizedPareto) = (d.ε, d.σ, d.ξ)
partype(::GeneralizedPareto{T}) where {T} = T

########################
### DENSITY FUNCTION ###
function pdf(d::GeneralizedPareto{T}, x::Real) where {T<:Real}
    p = -T(Inf)
    if x >= d.ε
        z = (x - d.ε) / d.σ
        if d.ξ > 0 || (d.ξ < 0 && x < maximum(d))
            p = 1/d.σ * (1 + d.ξ * z)^(-1 - 1/d.ξ)
        end
    end
    return p
end

function logpdf(d::GeneralizedPareto{T}, x::Real) where {T<:Real}
    (ε, σ, ξ) = params(d)

    # The logpdf is log(0) outside the support range.
    p = -T(Inf)

    if x >= ε
        z = (x - ε) / σ
        if abs(ξ) < eps()
            p = -z - log(σ)
        elseif ξ > 0 || (ξ < 0 && x < maximum(d))
            p = (-1 - 1 / ξ) * log1p(z * ξ) - log(σ)
        end
    end

    return p
end

#########################
### SURVIVAL FUNCTION ###
function logccdf(d::GeneralizedPareto, x::Real)
    ε, σ, ξ = params(d)
    z = max((x - ε) / σ, 0) # z(x) = z(ε) = 0 if x < ε (lower bound)
    return if abs(ξ) < eps(one(ξ)) # ξ == 0
        -z
    elseif ξ < 0
        # y(x) = y(ε - σ / ξ) = -1 if x > ε - σ / ξ (upper bound)
        -log1p(max(z * ξ, -1)) / ξ
    else
        -log1p(z * ξ) / ξ
    end
end
ccdf(d::GeneralizedPareto, x::Real) = exp(logccdf(d, x))

cdf(d::GeneralizedPareto, x::Real) = -expm1(logccdf(d, x))
logcdf(d::GeneralizedPareto, x::Real) = log1mexp(logccdf(d, x))

#################
### QUANTILES ###
function quantile(d::GeneralizedPareto{T}, p::Real) where T<:Real
    (ε, σ, ξ) = params(d)

    if p == 0
        z = zero(T)
    elseif p == 1
        z = ξ < 0 ? -1 / ξ : T(Inf)
    elseif 0 < p < 1
        if abs(ξ) < eps()
            z = -log1p(-p)
        else
            z = expm1(-ξ * log1p(-p)) / ξ
        end
    else
      z = T(NaN)
    end

    return ε + σ * z
end


################
### SAMPLING ###
function xval(d::GeneralizedPareto, u::Real)
    if abs(d.ξ) < eps()
        rd = -log(u)
    else
        rd = expm1(-d.ξ * log(u)) / d.ξ
    end
    return d.ε + d.σ * rd
end
rand(rng::AbstractRNG, d::GeneralizedPareto{T}) where {T<:Real} = xval(d, Random.rand(rng,float(T)))
function rand(rng::AbstractRNG, d::GeneralizedPareto{T}, n::Int) where {T<:Real}
    return map(Base.Fix1(xval, d), Random.rand(rng,float(T),n))
end
function rand!(rng::AbstractRNG, d::GeneralizedPareto, U::AbstractArray{<:Real})
    Random.rand!(rng, U)
    map!(Base.Fix1(xval, d), U, U)
    return U
end

###############
### FITTING ###
"""
Fit ParetoIV with for an array (Vector) of candidate minimum values εs
"""
function fit(::Type{GeneralizedPareto}, x::Array{T}, εs::Array{T}) where {T<:Real}
    xs = sort(x)

    σhat = eps()
    ξhat = eps()
    εhat = 1e-24
    D = Inf
    n = 0
    #/ For each possible xmin in xmins;
    #  - compute the max.-likelihood estimate of the power law exponent γ
    #  - compute the Kolmogorov-Smirnov distance
    #  - extract the xmin for which the MLE γ gives the smallest KS distance
    k = 0
    for i in eachindex(εs)
        ε = εs[i]
        n = count(xs .> ε)
        #~ Filter data
        _idx = searchsortedfirst(xs, ε)
        _x = xs[_idx:end]
        (length(_x) < 3) && (continue)
        _P = fit(GeneralizedPareto, _x, ε)
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        #  note: within the function data is filtered, so no need to do it here
        Dhat = KolmogorovSmirnov(_P, x)
        #~ If smaller than the current best, update
        if Dhat < D
            σhat = _P.σ
            ξhat = _P.ξ
            εhat = _P.ε
            D = Dhat
        end
        k += 1
    end
    return GeneralizedPareto(εhat, σhat, ξhat)
end

"""
    fit(d::GeneralizedPareto{T}, x::Array{T}; ε=0.0)

Fit a generalized Pareto distribution to data x>=ε [i.e., ε=xmin]
Uses method of moments as the initial guess for the parameters σ and ξ
Note that the mean vanishes when ξ<1, and so for real data fits with ξ>1 can be considered moot
"""
function fit(::Type{GeneralizedPareto}, x::Array{T}, ε::Float64) where {T<:Real}
    function negloglikelihood(x, θ)
        logσ, logξ = θ
        σ = exp(logσ)
        ξ = exp(logξ)
        d = GeneralizedPareto(ε, σ, ξ)
        #~ Return negative log-likelihood
        return -sum(logpdf.(d, x))
    end
    #~ Use method of moments for initial estimate
    Ex, Vx = StatsBase.mean(x), StatsBase.var(x)
    ξinit = (1 - (Ex - ε)^2 / Vx) / 2
    σinit = max(eps(), (Ex - ε)*(1-ξinit))
    θinit = [log(σinit), log(ξinit)]

    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, x),
        [log10(1e-10), log10(1e-8)],
        [log10(1e2), log10(3.0)],
        θinit,
        Fminbox(LBFGS()),
        Optim.Options(g_tol = 1e-3),
        autodiff=:forward
    )
    if Optim.converged(optimres)
        σhat, ξhat = optimres.minimizer
        return GeneralizedPareto(ε, exp(σhat), exp(ξhat))
    end
    throw(ErrorException("Optimizer not converged"))
    # println("Optimizer not converged, returning initial guesses [method of moments]")
end

########################
### HELPER FUNCTIONS ###
function KolmogorovSmirnov(P::GeneralizedPareto, data::Array{T}) where {T<:Real}
    #~ We care only about data within the functions domain, so first filter
    x = filter(z -> z >= P.ε, data)
    #~ Sort, if not already
    (!issorted(x)) && (sort!(x))
    Fv = _ecdf(x, x).F            # Values of empirical CDF
    Ftv = 1.0 .- ccdf.(P, x)      # Values of survival function
    KS = abs.(Fv .- Ftv)
    return maximum(KS)
end
