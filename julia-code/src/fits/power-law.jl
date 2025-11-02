#= Module for fitting Pareto-type distributions (power-laws) to data =#
#/ Start module
module Powerlaw

#/ Packages
using Distributions
using Random
using StatsBase

#/ Local modules
import Meris.ParetoLike as ParetoLike

#################
### FUNCTIONS ###
"""
Fit a power law on the heavy-tail of the data above some xmin
Uses methods from Clauset et al. (2009)
"""
function fitpowerlaw(x::Array{T}) where T<:Real
    xs = sort(x)
    γ = nothing
    xmin = 0.
    D = Inf
    n = 0
    for i in eachindex(xs)
        #~ Compute estimate of power-law exponent γ
        n = count(xs .> xs[i])
        (n < 50) && (break)        # If less than 50 samples >xmin, break
        _x = xs[i+1:end]
        _xmin = xs[i]
        S = sum(log.(_x / _xmin))
        γhat = 1 + n / S
        #~ Compute Kolmogorov-Smirnov distance
        Fv = _ecdf(_x, _x).F                    # Values of empirical CDF
        Ft = ParetoLike.Paretocdf(γhat, xmin=_xmin)  # Function of theoretical CDF
        Ftv = Ft.(_x)
        Z = sqrt.(Ftv .* (1 .- Ftv))
        distances = abs.(Fv .- Ftv) ./ Z
        Dhat = maximum(distances)
        if Dhat < D
            xmin = _xmin
            D = Dhat
            γ = γhat
        end
    end
    return (; γ=γ, xmin=xmin, KS=D)
end


"""
Bootstrap to compute p-values
"""
function bootstrap(
    x::Array{T}, γ::Float64, xmin::Float64, KS::Float64;
    rng=Random.Xoshiro(42),
    ε=0.1
    ) where T<:Real
	  nsynth = ceil(Int, 1 / (4*ε*ε))   #~ No. of synthetic datasets to generate
    #/ Compute the no. of samples in the tail
    n = length(x)
    xs = sort(x)
    xnonpowerlaw = xs[xs .<= xmin]
    ntail = n - searchsortedfirst(xs, xmin)
    ptail = ntail / n
    #/ Generate synthetic datasets and compute test statistics (Kolmogorov-Smirnov distance)
    pcounts = 0
    for m in 1:nsynth
        #~ Generate
        xsynth = similar(x, Float64)
        for k in 1:n
            utail, u = rand(rng, 2)
            if utail < ptail
                xsynth[k] = _samplepareto(u, γ, xmin)
            else
                xsynth[k] = StatsBase.sample(xnonpowerlaw)
            end
        end
        paretofit = fitpowerlaw(xsynth)
        if paretofit.KS > KS
            pcounts += 1
        end
    end
    return pcounts / nsynth
end


########################
### HELPER FUNCTIONS ###
########################
### DISTANCE FUNCTIONS
"""
Compute the Kolmogorov-Smirnov distance between the current best empirical CDF and theoretical CDF
"""
function KolmogorovSmirnov()
	  
end

### SAMPLERS
"""
Sample from Pareto distribution
"""
function samplepareto(
    nsamples::Int;
    γ::Float64 = 2.0,
    xmin::Float64 = 1e-4,
    rng = Random.Xoshiro(42*nsamples)
    )
	  u = rand(rng, nsamples)
    return _samplepareto.(u, γ, xmin)
end

function _samplepareto(u, γ, xmin)
	  return xmin .* (1 .- u).^(-1 / (γ .- 1.0))
end

"""
Sample from the Burr distribution

# Notes:
- when `c=1`, the Burr distribution is the Lomax distribution
- when `α=1`, the Burr distribution is a log-logistic distribution
"""
function sampleburr(
    nsamples::Int;
    c::Float64=1.0,
    α::Float64=1.0,
    λ::Float64=1.0,
    rng = Random.Xoshiro(42*nsamples)
    )
    u = rand(rng, nsamples) .^ (1 / α)
    return λ .* (1 ./ u .- 1) .^ (1 / c)    
end

"""
Sample from the Lomax distribution
"""
function samplelomax(
    nsamples::Int;
    α::Float64=1.0,
    λ::Float64=1.0,
    rng=Random.Xoshiro(42*nsamples)
    )
	  return sampleburr(nsamples; c=1.0, α=α, λ=λ, rng=rng)
end

"""
Compute empirical CDF at points t where F[t] = (no. elements ≤ t) / n
"""
function _ecdf(x::Array{T}, t::Array{T}; sorted=false) where T<:Real
    (!sorted) && (xs = sort(x))
    n = length(x)
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
    n = length(x)
    edges = range(xs[begin], xs[end], length=t) |> collect
    return _ecdf(x, edges, sorted=true)
end


"""
Compute empirical CDF at data points themselves
"""
function _ecdf(x::Array{T}) where T<:Real
    return _ecdf(x, sort(x), sorted=false)
end

######################
### TEST FUNCTIONS ###
"""
Do the same analysis on a bunch of synthetic data set, and compute p-value statistics
"""
function evaluate_powerlawfit(nsamples, nsims; γ=2.5, xmin=1.0)
    ps = Float64[]
    for i in 1:nsims
        x = samplepareto(nsamples; γ=γ, xmin=xmin, rng=Random.Xoshiro(42 + i))
        fit = fitpowerlaw(x)
        p = bootstrap(x, fit.γ, fit.xmin, fit.KS; rng=Random.Xoshiro(12 + i))
        push!(ps, p)
    end
    return ps
end

end # module Powerlaw
#/ End module
