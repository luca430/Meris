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
Fit ParetoIV with for an array (Vector) of candidate minimum values εs
"""
function fit(::Type{ParetoIV}, x::Array{T}, εs::Vector{T}; weighted=false) where {T<:Real}
    (!issorted(x)) && (x = sort(x))

    αhat = eps()
    βhat = eps()
    θhat = eps()
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
        #~ Filter data
        _idx = searchsortedfirst(x, ε)
        _x = x[_idx:end]
        n = length(_x)
        (n < 32) && (continue)
        _P = fit(ParetoIV, _x, ε)
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        #  note: within the function data is filtered, so no need to do it here
        Dhat = KolmogorovSmirnov(_P, x; weighted=weighted)
        #~ If smaller than the current best, update
        if Dhat < D
            αhat = _P.α
            βhat = _P.β
            θhat = _P.θ
            εhat = _P.ε
            D = Dhat
        end
        k += 1
    end
    return ParetoIV(αhat, βhat, θhat, εhat)
end

"""
Fit ParetoIV with fixed inequality paremeter β for an array (Vector) of candidate
minimum values εs. Fixing β allows to specify fitting of other functional forms.
For example, when β=1 then the ParetoIV is a Burr/Lomax distribution.
"""
function fit(::Type{ParetoIV}, x::Array{T}, β::T, εs::Array{T}; weighted=false) where {T<:Real}
    (!issorted(x)) && (x = sort(x))

    αhat = eps()
    θhat = eps()
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
        #~ Filter data
        _idx = searchsortedfirst(x, ε)
        _x = x[_idx:end]
        n = length(_x)
        (n < 32) && (break)
        _P = fit(ParetoIV, _x, β, ε)
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        #  note: within the function data is filtered, so no need to do it here
        Dhat = KolmogorovSmirnov(_P, x; weighted=weighted)
        #~ If smaller than the current best, update
        if Dhat < D
            αhat = _P.α
            θhat = _P.θ
            εhat = _P.ε
            D = Dhat
        end
        k += 1
    end
    return ParetoIV(αhat, β, θhat, εhat)
end

"""
Fit ParetoIV with a fixed minimum value ε
"""
function fit(::Type{ParetoIV}, x::Array{T}, ε::T) where {T<:Real}
    #/ Negative log-likelihood
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
    θinit = iqr(x .- ε)
    #~ Simple ParetoI MLE as initial estimate
    _x = x[x.>ε]
    S = sum(log.(_x / ε))
    αinit = max(0.1, length(_x) / S)
    #~ Guess of β, typically β∈[0.5,1.0]
    βinit = .99
    params = [log(αinit), log(βinit), log(θinit)]
    #~ Optimize
    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, x),
        [log(1e-3), log(1e-3), log(1e-8)],
        [log(10.0), log(1.0), log(maximum(x))],
        params,
        Fminbox(LBFGS()),
        autodiff=:forward
    )
    if Optim.converged(optimres)
        αhat, βhat, θhat = optimres.minimizer
        return ParetoIV(exp(αhat), exp(βhat), exp(θhat), ε)
    end
    #~ Throw an error here as the Optimizer did not converge
    throw(ErrorException("Optimizer not converged"))
    #~ When an error is undesired, it can be set as a warning as well
    # @warn("Optimizer not converged, returning initial guesses")
    # return ParetoIV(αinit, βinit, θinit, ε)
end

"""
Fit ParetoIV with a fixed inequality parameter β and fixed minimum value ε
"""
function fit(::Type{ParetoIV}, x::Array{T}, β::T, ε::T) where {T<:Real}
    #/ Negative log-likelihood
    function negloglikelihood(x, params)
        logα, logθ = params
	      α = exp(logα)
        θ = exp(logθ)
        d = ParetoIV(α, β, θ, ε)
        return -sum(logpdf.(d, x))
    end

    #/ Init. estimates
    #~ Quantile estimator for θ
    θinit = iqr(x .- ε)
    #~ Simple ParetoI MLE as initial estimate
    _x = x[x.>ε]
    S = sum(log.(_x / ε))
    αinit = max(0.1, length(_x) / S)
    params = [log(αinit), log(θinit)]
    #~ Optimize
    optimres = Optim.optimize(
        Base.Fix1(negloglikelihood, _x),
        [log(1e-3), log(1e-6)],
        [log(10.0), log(maximum(x))],
        params,
        Fminbox(LBFGS()),
        Optim.Options(g_tol = 1e-3),
        autodiff=:forward
    )
    if Optim.converged(optimres)
        αhat, θhat = optimres.minimizer
        return ParetoIV(exp(αhat), β, exp(θhat), ε)
    end
    #~ Throw an error here as the Optimizer did not converge
    throw(ErrorException("Optimizer not converged"))
    #~ When an error is undesired, it can be set as a warning as well
    # @warn("Optimizer not converged, returning initial guesses")
    # return ParetoIV(αinit, βinit, θinit, ε)
end

"""
Compute p-value that determines whether to reject the ParetoI as a candidate
see, [Clauset et al. (2009), Power-law distribution in empirical data]
"""
function computepvalue(
    P::ParetoIV, x::Array{T}, εs::Array{T};
    nsynth=625, weighted=false, rng=Random.Xoshiro(42)
    ) where {T<:Real}
    #~ Compute prob. to augment synthetic data [see Clauset et al. (2009), Section 4.1]
    xhead = filter(z -> z < P.ε, x)
    nhead = length(xhead) / length(x)
    k = length(x)
    #~ Compute Kolmogorov-Smirnov distance in data
    KSDATA = KolmogorovSmirnov(P, x; weighted=weighted)
    kscount = 0

    #/ Generate synthetic datasets
    ns = 0
    while ns < nsynth
        ns += 1
        #~ Sample synthetic dataset
        synthx = rand(rng, P, k)
        #~ augment synthetic data [see Clauset et al. (2009), Section 4.1]
        u = Base.rand(rng, length(x))
        for i in eachindex(synthx)
            (u[i] < nhead) && (synthx[i] = StatsBase.sample(rng, xhead))
        end        
        #~ Choose admissible ε
        logsynthx = log.(synthx)
        logxmin, logxmax = extrema(logsynthx)
        εsynth = exp.(range(logxmin, logxmax, length(εs)))
        Psynthfit = fit(ParetoIV, synthx; εs=εsynth, weighted=weighted)
        #~ Compute Kolmogorov-Smirnov distance in synthetic data
        KSSYNTHETIC = KolmogorovSmirnov(Psynthfit, synthx; weighted=weighted)
        if KSSYNTHETIC > KSDATA
            kscount += 1
        end
    end
    #~ return p value
    return kscount / nsynth
end

"""
Compute p-value for fixed β in order to determine whether to reject the ParetoIV as a candidate
see, [Clauset et al. (2009), Power-law distribution in empirical data]
"""
function computepvalue(
    P::ParetoIV, x::Vector{T}, β::T, εs::Vector{T};
    nsynth=332, weighted=false
) where {T<:Real}
    #~ Compute prob. of data appearing in the "head" [the "non-tail"]
    xhead = filter(z -> z < P.ε, x)
    nhead = length(xhead) / length(x)
    k = length(x)
    KSDATA = KolmogorovSmirnov(P, x; weighted=weighted)
    #~ Allocate parallel thread counts and rngs
    kscount = Threads.Atomic{Int}(0)
    rngs = [Random.Xoshiro(42*t) for t in 1:Threads.nthreads()]

    Threads.@threads for _ in 1:nsynth
        #~ Generate synthetic data using the proper rng
        _rng = rngs[Threads.threadid()]
        synthx = rand(_rng, P, k)
        #~ augment synthetic data [see Clauset et al. (2009), Section 4.1]
        u = Base.rand(_rng, k)
        for i in eachindex(synthx)
            (u[i] < nhead) && (synthx[i] = StatsBase.sample(_rng, xhead))
        end
        #~ Choose admissible ε
        logsynthx = log.(synthx)
        logxmin, logxmax = extrema(logsynthx)
        εsynth = exp.(range(logxmin, logxmax, length(εs)))
        #~ Fit on synthetic data
        Psynthfit = fit(ParetoIV, synthx, β, εsynth; weighted=weighted)
        #~ Compute Kolmogorov-Smirnov distance in synthetic data
        KSSYNTHETIC = KolmogorovSmirnov(Psynthfit, synthx; weighted=weighted)
        #~ Compute Kolmogorov-Smirnov distance in synthetic data and, if larger, increment
        if KSSYNTHETIC > KSDATA
            Threads.atomic_add!(kscount, 1)
        end
    end
    #~ return p value
    return kscount[] / nsynth
end

########################
### HELPER FUNCTIONS ###
function KolmogorovSmirnov(P::ParetoIV, data::Array{T}; weighted=false) where {T<:Real}
    #~ We care only about data within the functions domain, so first filter
    x = filter(z -> z >= P.ε, data)
    #~ Sort, if not already
    (!issorted(x)) && (sort!(x))
    Fv = _ecdf(x, x).F            # Values of empirical CDF
    Ftv = 1.0 .- ccdf.(P, x)      # Values of survival function
    if weighted
        Z = sqrt.(Ftv .* (1 .- Ftv))    # Weight
        KS = abs.(Fv .- Ftv) ./ Z       # Weighted KS distance
        return maximum(KS)
    end
    KS = abs.(Fv .- Ftv)
    return maximum(KS)
end

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
