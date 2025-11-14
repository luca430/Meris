#= Module for fitting Pareto-type distributions (power-laws) to data =#
#/ Start module
module Powerlaw

#/ Packages
using Distributions
using Random
using StatsBase
using SpecialFunctions
using RootSolvers

using FHist

#/ Local modules
using Meris

#################
### FUNCTIONS ###
"""
Fit a power law on the heavy-tail of the data above some xmin
Uses methods from Clauset et al. (2009)
"""
function fitPareto(x::Array{T}; xmins=nothing) where T<:Real
    xs = sort(x)
    xmins = isnothing(xmins) ? unique(xs) : xmins
    γ = nothing
    xmin = 0.
    D = Inf
    n = 0
    #/ For each possible xmin in xmins;
    #  - compute the max.-likelihood estimate of the power law exponent γ
    #  - compute the Kolmogorov-Smirnov distance
    #  - extract the xmin for which the MLE γ gives the smallest KS distance
    for i in eachindex(xmins)
        _xmin = xmins[i]
        n = count(xs .>= _xmin)
        (n < 50) && (break)        # If less than 50 samples >xmin, break
        #~ Filter data
        _idx = searchsortedfirst(xs, _xmin)
        _x = xs[_idx:end]        
        #~ Compute MLE of power-law exponent γ
        S = sum(log.(_x / _xmin))
        γhat = 1 + n / S
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        Fv = _ecdf(_x, _x, sorted=true).F                  # Values of empirical CDF
        Ft = Meris.ParetoLike.Paretocdf(γhat, xmin=_xmin)  # Function of theoretical CDF
        Ftv = Ft.(_x)                                      # Values of the theoretical CDF
        Z = sqrt.(Ftv .* (1 .- Ftv))                       # Weighted KS distance

        #~ @TODO Anderson-Darling test statistic
        #  Clauset et al. note that this test statistic may be 'too' conservative, especially
        #  where there are not 'enough' samples in the tail.
        # s = 1:n
        # S = sum((2 .* s .- 1) ./ n .* (log.(Ftv) .+ log.(1 .- reverse(Ftv))))
        # Dhat = -(n+S)        
        
        distances = abs.(Fv .- Ftv) ./ Z
        Dhat = maximum(distances)
        #~ If smaller than the current best, update
        if Dhat < D
            xmin = _xmin
            D = Dhat
            γ = γhat
        end
    end
    return (; γ=γ, xmin=xmin, KS=D)
end

"""
Fit a power law on the heavy-tail of the data above some xmin
Uses methods from Clauset et al. (2009)
"""
function fitGeneralizedPareto(x::Array{T}; xmins=nothing) where T<:Real
    xs = sort(x)
    xmins = isnothing(xmins) ? unique(xs) : xmins
    σ = nothing
    ξ = nothing
    xmin = 0.
    D = Inf
    n = 0
    #/ For each possible xmin in xmins;
    #  - compute the max.-likelihood estimate of the power law exponent γ
    #  - compute the Kolmogorov-Smirnov distance
    #  - extract the xmin for which the MLE γ gives the smallest KS distance
    for i in eachindex(xmins)
        _xmin = xmins[i]
        n = count(xs .>= _xmin)
        (n < 50) && (break)        # If less than 50 samples >xmin, break
        #~ Filter data
        _idx = searchsortedfirst(xs, _xmin)
        _x = xs[_idx:end]        
        #~ Compute MLE of generalized Pareto distribution with μ=xmin
        gpdparams = Meris.MLEFit.fitgeneralizedPareto(_x; μ=_xmin)
        #  note: assumes convergence of Optim, otherwise good results cannot be guaranteed
        σhat = gpdparams.σ
        ξhat = gpdparams.ξ
        
        #~ Compute Kolmogorov-Smirnov distance as the test statistic
        Fv = _ecdf(_x, _x, sorted=true).F                  # Values of empirical CDF
        Ft = Meris.ParetoLike.generalizedParetocdf(σhat, ξhat, xmin=_xmin)
        Ftv = Ft.(_x)                                      # Values of the theoretical CDF
        Z = sqrt.(Ftv .* (1 .- Ftv))                       # Weighted KS distance
        
        distances = abs.(Fv .- Ftv) ./ Z
        Dhat = maximum(distances)
        
        #~ If smaller than the current best, update
        if Dhat < D
            xmin = _xmin
            D = Dhat
            σ = σhat
            ξ = ξhat
        end
    end
    return (; σ=σ, ξ=ξ, xmin=xmin, KS=D)
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

"""
    log_variance(x::Array{T})

Compute the log-variance using methods from Lee & Kim (2018) to estimate ξ = 1 / α
[for details, see: https://doi.org/10.1080/03610926.2018.1441418]
"""
function log_variances(x::Array{T}; nmin::Int=32, sorted=false) where T<:Real
    n = length(x)
	  xs = sorted ? x : reverse(sort(x))          #~ Sort descending
    ξ = zeros(Float64, n - nmin)

    for i in 1:(n-nmin)
        #~ Gather exceedences xs > xs[i]
        idx = searchsortedfirst(xs, xs[i], lt = >)
        (idx < nmin) && (continue)
        #~ Compute log-exceedences
        logx = log.(xs[begin:idx-1] .- xs[i])
        #~ Estimate ξ
        #  note: The inversetrigamma is estimated here with a few Newton steps
        #        as SpecialFunctions does not have this implemented directly.
        s = var(logx)
        ξ[i] = 1 / invtrigamma(s - trigamma(1))
    end
    return ξ
end


"""
Bootstrap to compute p-values
"""
function bootstrapPareto(
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
            #~ With prob. ptail, generate Pareto samples,
            #  otherwise take a sample from the data with x<xmin
            if utail < ptail
                xsynth[k] = _samplePareto(u, γ, xmin)
            else
                xsynth[k] = StatsBase.sample(xnonpowerlaw)
            end
        end
        paretofit = fitPareto(xsynth)
        if paretofit.KS > KS
            pcounts += 1
        end
    end
    return pcounts / nsynth
end


########################
### HELPER FUNCTIONS ###
########################
### SPECIAL FUNCTIONS
"""
    invtrigamma(x::Float64)

Gives inverse trigamma, x = (ψ')^{-1}(y), where y = ψ'(x)
As the trigamma function is strictly decreasing for x>0, it is invertable and we can use a
simple/naive root-finding method, such as Newton-Raphson, to get a good estimate of its inverse.
Typically needs only a handful of iterations for the required tolerance.
"""
function invtrigamma(y; tol=1e-12, maxiter=32)
    #~ Make a good initial guess [from Taylor expansion]
    x = 1 ./ y .- 0.5
    #~ Allocate
    iter = 0
    τ = Inf
    #~ Iterate simple Newton-Raphson until convergence or max iterations
    while τ > tol && iter < maxiter
        dx = (trigamma.(x) .- y) ./ polygamma.(2, x)
        τ = abs(dx)
        x -= dx
        iter += 1
    end
    return x
end

### SAMPLERS
"""
Sample from Pareto distribution

A Pareto distribution has survival function S(x) = c x⁻ᵅ
Uses inverse transform sampling as the inverse has a simple closed form
"""
function samplePareto(
    nsamples::Int;
    α::Float64 = 2.0,
    xmin::Float64 = 1e-4,
    rng = Random.Xoshiro(42*nsamples)
    )
	  u = rand(rng, nsamples)
    return _samplePareto.(u, α, xmin)
end

function _samplePareto(u, α, xmin)
	  return xmin .* (1 .- u).^(-1 / (α .- 1.0))
end

"""
Sample from tempered Pareto distribution

A tempered Pareto distribution has survival function S(x) = c x⁻ᵅ e⁻ᵝˣ
"""
function sampletemperedPareto(
    nsamples::Int;
    α::Float64 = 2.0,
    β::Float64 = 0.0,
    xmin::Float64 = 1e0,
    rng = Random.Xoshiro(42*nsamples)
)
    (iszero(β)) && (return samplePareto(nsamples, α=α, xmin=xmin, rng=rng))
    u = rand(rng, nsamples)
    return map(r -> _sampletemperedPareto(r, α, β, xmin), u)
end

function _sampletemperedPareto(u, α, β, xmin)
    z = log(xmin^α * exp(β*xmin) / (1 - u))
    #~ Define (function, deriviative) for NewtonsMethod
    #! note: we do y = log(x) to avoid Newton's method to have x < 0 [see below]
    fdf(x) = (α*x + β*exp(x) - z, α + β*exp(x))
    #~ Solve with good initial guess
    #! note: as x>0, we need the initial guess to be "close", as otherwise the method may
    #        overshoot into x<0 territory, so the guess must depend on `u`.
    guess = log(max(xmin, z / β))
    sol = RootSolvers.find_zero(fdf, RootSolvers.NewtonsMethod{Float64}(guess))
    #~ As we solved for y = log(x), exponentiate
    return exp(sol.root)
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
    u = rand(rng, nsamples)
    return λ .* ((1 .- u).^(-1/α) .- 1) .^ (1 / c)
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

#
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

######################
### TEST FUNCTIONS ###
"""
Do the same analysis on a bunch of synthetic data set, and compute p-value statistics
"""
function evaluate_powerlawfit(nsamples, nsims; γ=2.5, xmin=1.0)
    ps = Float64[]
    for i in 1:nsims
        x = samplePareto(nsamples; γ=γ, xmin=xmin, rng=Random.Xoshiro(42 + i))
        fit = fitpowerlaw(x)
        p = bootstrap(x, fit.γ, fit.xmin, fit.KS; rng=Random.Xoshiro(12 + i))
        push!(ps, p)
    end
    return ps
end

function checkburr(nsamples; c::Float64=1.5, α::Float64=0.5, λ::Float64=2.0)
    xmin = -5
    xmax = 5
	  #~ Sample from Burr
    xs = sampleburr(nsamples; c=c, α=α, λ=λ)
    bins = exp10.(range(xmin,xmax,32))
    fh = FHist.Hist1D(xs; binedges=bins) |> FHist.normalize

    xburr = exp10.(range(xmin,xmax,256)) |> collect
    F = Meris.ParetoLike.Burrcdf(c, α, λ)
    yburr = Meris.ParetoLike.Burrpdf(xburr, c, α, λ) .* (1 - F(xmax))
    
    return (; x=bincenters(fh), y=fh.bincounts), (; x=xburr, y=yburr)
end

end # module Powerlaw
#/ End module
