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
    βinit = max(2/maximum(x), 1 / (Ex + ε))
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

function fit(::Type{TemperedPareto}, x::Array{T}; εs=nothing, weighted=false) where {T<:Real}
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
        #~ Filter data
        _idx = searchsortedfirst(xs, ε)
        _x = xs[_idx:end]
        n = length(_x)
        (n < 64) && (break)    #~ `break` if not enough samples remain
        _P = fit(TemperedPareto, _x, ε)
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        #  note: within the function data is filtered, so no need to do it here
        Dhat = KolmogorovSmirnov(_P, x; weighted=weighted)
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

"""
Compute p-value that determines whether to reject the tempered Pareto as a candidate
see, [Clauset et al. (2009), Power-law distribution in empirical data]
"""
function computepvalue(
    P::TemperedPareto, x::Array{T}, εs::Array{T};
    nsynth = 1000, weighted=false, rng=Random.Xoshiro(42)
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
        Psynthfit = fit(ParetoI, synthx; εs=εsynth, weighted=weighted)
        #~ Compute Kolmogorov-Smirnov distance in synthetic data
        KSSYNTHETIC = KolmogorovSmirnov(Psynthfit, synthx; weighted=weighted)
        if KSSYNTHETIC > KSDATA
            kscount += 1
        end
    end
    #~ return p value
    return kscount / nsynth
end


########################
### HELPER FUNCTIONS ###
function KolmogorovSmirnov(P::TemperedPareto, data::Array{T}; weighted=false) where {T<:Real}
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
