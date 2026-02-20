###############
### STRUCTS ###

struct ParetoI{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    ε::T

    ParetoI{T}(α::T, ε::T) where {T<:Real} = new{T}(α, ε)
end

function ParetoI(α::T, ε::T; check_args::Bool = true) where {T<:Real}
	  @check_args ParetoI (α, α>zero(α)) (ε, ε>zero(ε))
    return ParetoI{T}(α,ε)
end

####################
### CONSTRUCTORS ###

Pareto(α::Real,ε::Real; check_args::Bool=true) = ParetoI(promote(α,ε)...; check_args=check_args)
ParetoI(α::Real,ε::Real; check_args::Bool=true) = ParetoI(promote(α,ε)...; check_args=check_args)

##################
### STATISTICS ###

mean(d::ParetoI) = if d.α <= 1 return Inf else return d.α*d.ε / (d.α - 1) end
scale(d::ParetoI) = d.ε
shape(d::ParetoI) = d.α
params(d::ParetoI) = (d.α, d.ε)

########################
### DENSITY FUNCTION ###

function pdf(d::ParetoI, x::Real)
	  if x < d.ε
        return 0.0
    end
    return d.α*d.ε^d.α * x^(-(d.α+1))
end

function logpdf(d::ParetoI, x::Real)
    (x < d.ε) && (return -Inf)
    return d.α * log(d.α * d.ε) - (d.α + 1) * log(x)
end

#########################
### SURVIVAL FUNCTION ###

function ccdf(d::ParetoI, x::Real)
	  if x < d.ε
        return 1.0
    end
    return (d.ε / x)^d.α
end

################
### SAMPLING ###

function xval(d::ParetoI, u::Real)
    return d.ε / u^(1/d.α)
end
rand(rng::AbstractRNG, d::ParetoI{T}) where {T<:Real} = xval(d, Random.rand(rng,float(T)))
function rand(rng::AbstractRNG, d::ParetoI{T}, n::Int) where {T<:Real}
    return map(Base.Fix1(xval, d), Random.rand(rng,float(T),n))
end
function rand!(rng::AbstractRNG, d::ParetoI, U::AbstractArray{<:Real})
    Random.rand!(rng, U)
    map!(Base.Fix1(xval, d), U, U)
    return U
end

###############
### FITTING ###

function fit(::Type{ParetoI}, x::Array{T}, ε::Float64) where {T<:Real}
    (!issorted(x)) && (x = sort(x))
    n = count(x .>= ε)
    (n < 256) && (@warn("Very little data in the tail with ε=$(ε) [only $(n) data points]"))
    #~ Filter data
    idx = searchsortedfirst(x, ε)
    xfit = x[idx:end]
    #/ Compute maximum-likelihood estimate of the ParetoI exponent α
    S = sum(log.(xfit / ε))
    αhat = n / S
    return ParetoI(αhat, ε)
end

function fit(::Type{ParetoI}, x::Array{T}; εs=nothing, weighted=false) where {T<:Real}
    xs = sort(x)
    εs = isnothing(εs) ? unique(xs) : εs
    
    αhat = 1.0
    εhat = 1e-24
    D = Inf
    n = 0
    #/ For each possible xmin in xmins;
    #  - compute the max.-likelihood estimate of the power law exponent γ
    #  - compute the Kolmogorov-Smirnov distance
    #  - extract the ε (xmin) for which the MLE γ gives the smallest KS distance
    for i in eachindex(εs)
        ε = εs[i]
        #~ Filter data
        _idx = searchsortedfirst(xs, ε)
        _x = xs[_idx:end]
        n = length(_x)
        (n < 128) && (break)        # If less than 256 samples >xmin, break
        #~ Estimate α
        S = sum(log.(_x / ε))
        α = n / S
        _P = ParetoI(α, ε)
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        #  note: within the function data is filtered, so no need to do it here
        Dhat = KolmogorovSmirnov(_P, x; weighted=weighted)
        #~ If smaller than the current best, update
        if Dhat < D
            αhat = _P.α
            εhat = _P.ε
            D = Dhat
        end
    end
    return ParetoI(αhat, εhat)
end

"""
Compute p-value that determines whether to reject the generalized Pareto as a candidate
see, [Clauset et al. (2009), Power-law distribution in empirical data]
"""
function computepvalue(
    P::ParetoI, x::Array{T}, εs::Array{T};
    nsynth=500, weighted=false, rng=Random.Xoshiro(42)
    ) where {T<:Real}
    #~ Compute prob. to augment synthetic data [see Clauset et al. (2009), Section 4.1]
    xhead = filter(z -> z < P.ε, x)
    nhead = length(xhead) / length(x)
    #~ Filter data
    # xs = sort(x[x .>= P.ε])
    k = length(x)
    #~ Compute Kolmogorov-Smirnov distance in data
    KSDATA = KolmogorovSmirnov(P, x[x .> P.ε]; weighted=weighted)
    kscount = 0

    #/ Generate synthetic datasets
    ns = 0
    while ns < nsynth
        ns += 1
        #~ Sample synthetic dataset
        synthx = sort(rand(rng, P, k))
        #~ augment synthetic data [see Clauset et al. (2009), Section 4.1]
        u = Base.rand(rng, length(x))
        for i in eachindex(x)
            (u[i] < nhead) && (synthx[i] = StatsBase.sample(rng, xhead))
        end        
        #~ Choose admissible ε
        _logx = log.(synthx)
        # εsynth = exp.(range(_logx[begin], _logx[end], length(εs)))
        Psynthfit = fit(ParetoI, synthx, P.ε)
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
function KolmogorovSmirnov(P::ParetoI, x::Array{T}; weighted=false) where {T<:Real}
    #~ We care only about data within the functions domain, so first filter
    filter!(z -> z >= P.ε, x)
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
Compute empirical CDF at points t where F[t] = (no. elements ≤ t) / n
"""
function _ecdf(xs::Array{T}, t::Array{T}) where T<:Real
    #~ ensure `t` is sorted
    (!issorted(t))  && (sort!(t))
    (!issorted(xs)) && (xs = sort(xs))
    n = length(xs)
    F = similar(t, Float64)
    k = 1
    for i in eachindex(t)
        #~ Move k until xs[k] > edges[k]
        while k ≤ n && xs[k] ≤ t[i]
            k += 1
        end
        F[i] = (k-1) / n
    end
    return (; F=F, t=t)
end

"""
Compute empirical CDF at equally distributed points t
"""
function _ecdf(xs::Array{T}, t::Int; sorted=false) where T<:Real
    (!sorted) && (xs = sort(xs))
    edges = range(xs[begin], xs[end], length=t) |> collect
    return _ecdf(xs, edges, sorted=true)
end
