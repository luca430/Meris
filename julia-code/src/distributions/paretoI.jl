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
    xs = sort(x)
    n = count(x .>= ε)
    (n < 256) && (@warn("Very little data in the tail with ε=$(ε) [only $(n) data points]"))
    #~ Filter data
    idx = searchsortedfirst(xs, ε)
    xfit = xs[idx:end]
    #/ Compute maximum-likelihood estimate of the ParetoI exponent α
    S = sum(log.(xfit / ε))
    αhat = n / S
    return ParetoI(αhat, ε)
end

function fit(::Type{ParetoI}, x::Array{T}; εs=nothing) where {T<:Real}
    xs = sort(x)
    εs = isnothing(εs) ? unique(xs) : εs
    
    αhat = 1.0
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
        (n < 256) && (break)        # If less than 256 samples >xmin, break
        S = sum(log.(_x / ε))
        α = n / S
        _P = ParetoI(α, ε)
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        Fv = _ecdf(_x, _x, sorted=true).F     # Values of empirical CDF
        Ftv = 1.0 .- ccdf.(_P, _x)
        Z = sqrt.(Ftv .* (1 .- Ftv))          # Weight
        distances = abs.(Fv .- Ftv) ./ Z      # Weighted KS distance
        Dhat = Base.maximum(distances)
        #~ If smaller than the current best, update
        if Dhat < D
            αhat = _P.α
            εhat = _P.ε
            D = Dhat
        end
    end
    @info "hm" D
    return ParetoI(αhat, εhat)
end

"""
Compute p-value that determines whether to reject the generalized Pareto as a candidate
see, [Clauset et al. (2009), Power-law distribution in empirical data]
"""
function computepvalue(P::ParetoI, x::Array{T}; nsynth=500) where {T<:Real}
    #~ Simulate `nsynth` synthetic datasets, and fit a generalized Pareto on each
    #  Then, compute the (weighted) KS statistic and compare with the value for the data
    rng = Random.Xoshiro(42*nsynth)
    
    #~ Compute KS statistic in data
    xs = sort(x[x .>= P.ε])
    k = length(xs)
    Fv = _ecdf(xs, xs, sorted=true).F     # Values of empirical CDF
    Ftv = 1.0 .- ccdf.(P, xs)
    # Z = sqrt.(Ftv .* (1 .- Ftv))          # Weight
    # distances = abs.(Fv .- Ftv) ./ Z      # Weighted KS distance
    distances = abs.(Fv .- Ftv)
    KSDATA = Base.maximum(distances)
    kscount = 0

    #/ Generate synthetic datasets
    for _ in 1:nsynth
        #~ Sample synthetic dataset
        _x = sort(rand(rng, P, k))
        #~ Fit a ParetoI with ε given
        S = sum(log.(_x / P.ε))
        α = k / S
        _P = ParetoI(α, P.ε)
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        Fv = _ecdf(_x, _x, sorted=true).F     # Values of empirical CDF
        Ftv = 1.0 .- ccdf.(_P, _x)
        # Z = sqrt.(Ftv .* (1 .- Ftv))          # Weight
        # distances = abs.(Fv .- Ftv) ./ Z      # Weighted KS distance
        distances = abs.(Fv .- Ftv)
        KSSYNTHETIC = Base.maximum(distances)
        if KSSYNTHETIC > KSDATA
            kscount += 1
        end
    end
    return kscount / nsynth
end

########################
### HELPER FUNCTIONS ###
"""
Compute empirical CDF at points t where F[t] = (no. elements ≤ t) / n
"""
function _ecdf(xs::Array{T}, t::Array{T}; sorted=false) where T<:Real
    #~ ensure `t` is sorted
    (!issorted(t)) && (sort!(t))
    (!sorted) && (xs = sort(xs))
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


"""
Compute empirical CDF at data points themselves
"""
function _ecdf(x::Array{T}) where T<:Real
    return _ecdf(x, sort(x), sorted=false)
end
