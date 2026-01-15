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

#########################
### SURVIVAL FUNCTION ###

function ccdf(d::ParetoI, x::Real)
	  if x <= d.ε
        return 1.0
    end
    return (d.ε / x)^d.α
end

################
### SAMPLING ###

function xval(d::ParetoI, u::Real)
    return d.ε / u^(1/(1-d.α))
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
    (n < 30) && (@warn("Very little data in the tail with ε=$(ε) [only $(n) data points]"))
    #~ Filter data
    idx = searchsortedfirst(xs, ε)
    xfit = xs[idx:end]
    #/ Compute maximum-likelihood estimate of the ParetoI exponent α
    S = sum(log.(xfit / ε))
    αhat = n / S
    return ParetoI(αhat, ε)
end

function fit(::Type{ParetoI}, x::Array{T}; εs=nothing) where T<:Real
    xs = sort(x)
    εs = isnothing(εs) ? unique(xs) : εs
    
    αhat = nothing
    εhat = nothing
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
        _P = fit(ParetoI, _x, ε)
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        Fv = _ecdf(_x, _x, sorted=true).F     # Values of empirical CDF
        Ftv = 1.0 .- ccdf.(_P, _x)
        Z = sqrt.(Ftv .* (1 .- Ftv))          # Weight
        distances = abs.(Fv .- Ftv) ./ Z      # Weighted KS distance
        Dhat = maximum(distances)
        #~ If smaller than the current best, update
        if Dhat < D
            αhat = _P.α
            εhat = _P.ε
            D = Dhat
        end
    end
    return ParetoI(αhat, εhat)
end

########################
### HELPER FUNCTIONS ###
"""
Compute empirical CDF at points t where F[t] = (no. elements ≤ t) / n
"""
function _ecdf(xs::Array{T}, t::Array{T}; sorted=false) where T<:Real
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
function _ecdf(x::Array{T}, t::Int; sorted=false) where T<:Real
    (!sorted) && (xs = sort(x))
    edges = range(xs[begin], xs[end], length=t) |> collect
    return _ecdf(x, edges, sorted=true)
end


"""
Compute empirical CDF at data points themselves
"""
function _ecdf(x::Array{T}) where T<:Real
    return _ecdf(x, sort(x), sorted=false)
end
