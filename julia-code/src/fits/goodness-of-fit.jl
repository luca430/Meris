#= Module to perform goodness-of-fit evaluations to check how well some data
   - is described by a particular function
   - is described by a particular distribution
   - is described by some distribution compared to another distribution
=#
#/ Start module
module GOF

#/ Packages
using Random
using FHist
using Statistics, StatsBase, LinearAlgebra

#################
### FUNCTIONS ###
"""
    estimatep(data, CDF, invCDF, test; nruns=10000, rng=Random.Xoshiro(42), params...)

Estimate p-value for a goodness-of-fit (GOF) test using Monte-Carlo simulations.
This function is intended specifically to perform a GOF test for distributions using one of the available methods.

# Arguments
- `data`: Sample data.
- `G`: Theoretical cumulative distribution function (CDF).
- `G_inv`: Inverse CDF (quantile function) used for sampling under the null hypothesis.
- `test`: A function implementing the test statistic
          (e.g., `chisquared`, `AndersonDarling`, `CramervonMises`, `KolmogorovSmirnov`).
- `nruns`: Number of Monte Carlo runs (default = 10_000).
- `rng`: Random number generator (default = `Random.Xoshiro(42)`).
- `params...`: Additional keyword arguments passed to the test function.

# Returns
Estimated p-value `p = mean(T_sim ≥ T)` where `T` is the observed statistic.
"""
function estimatep(data, CDF, invCDF, test; nmcsamples=10_000, rng=Random.Xoshiro(42), params...)
    n = length(data)
    #~ Compute the test statistics
    T = test(data, CDF; params...)
    #~ Generate samples from the inverse CDF [if available]
    samples = [invCDF.(rand(rng, Float64, n)) for _ in 1:nmcsamples]
    #~ Compute test statistics from the Monte-Carlo runs
    Tsim = [test(s, CDF; params...) for s in samples]
    #~ Compute the p-value
    p = mean(Tsim .>= T)
    return p
end

"""
    fit_scores(fit, model, x, y)

Estimate many statistic scores for a goodness-of-fit (GOF) test.
This function is intended to perform a GOF test for a generic model, evaluating most of the classical GOF estimators.

# Arguments
- `model`: Theoretical model in used for the fit. It must be in the form 'model(x, p)'.
- `x`: x data used to fit.
- `y`: y data used to fit.
- `p`: Array of model parameters

# Returns
Estimated p-value `p = mean(T_sim ≥ T)` where `T` is the observed statistic.
"""
function fit_scores(model, x, y, p)
    ŷ = model(x, p)
    n = length(y)
    k = length(p)

    # Residuals and SSE
    res = y .- ŷ
    SSE = sum(abs2, res)
    MSE = SSE / n
    RMSE = sqrt(MSE)

    # R^2 and adjusted R^2 (useful for regression-like fits)
    TSS = sum(abs2, y .- mean(y))
    R2 = 1 - SSE / TSS
    adjR2 = 1 - (1 - R2) * (n - 1) / (n - k - 1)

    # Estimate variance of residuals
    σ2_mle = SSE / n            # ML estimate (for AIC/BIC)
    σ2_unbiased = SSE / (n - k) # unbiased estimate (used for reduced chi^2)

    # Reduced chi-square (if you don't have per-point σ_i, use σ²_unbiased)
    χ2_reduced = σ2_unbiased == 0 ? Inf : (SSE / σ2_unbiased) / (n - k) # simplifies to 1, but keep form if you have explicit sigma_i

    # Log-likelihood under Gaussian iid errors (using ML variance = SSE/n)
    ll = -n/2 * (log(2π) + log(σ2_mle) + 1)
    AIC = 2k - 2ll
    BIC = k * log(n) - 2ll
    # AICc for small samples
    AICc = AIC + (2k*(k+1)) / (n - k - 1)

    return Dict(
        :n => n,
        :k => k,
        :SSE => SSE,
        :MSE => MSE,
        :RMSE => RMSE,
        :R2 => R2,
        :adjR2 => adjR2,
        :σ2_mle => σ2_mle,
        :σ2_unbiased => σ2_unbiased,
        :χ2_reduced => χ2_reduced,
        :loglikelihood => ll,
        :AIC => AIC,
        :AICc => AICc,
        :BIC => BIC,
        :residuals => res
    )
end


#######################
### TEST STATISTICS ###
"""
    chisquared(data, G; nbins=30)

Perform a χ² (chi-squared) goodness-of-fit test comparing `data` to distribution `G`.

# Arguments
- `data`:  Sample data.
- `CDF`:   Cumulative distribution function (CDF) of the theoretical distribution.
- `nbins`: Number of histogram bins [default = 30].

# Returns
The chi-squared test statistic.
"""
function chisquared(data, CDF; nbins::Int=30)
    #~ Create the histogram
    bmin, bmax = extrema(data)
    Δb = (bmax - bmin) / nbins
    fh = FHist.Hist1D(data, binedges=bmin:Δb:bmax, overflow=true)

    #~ Compute the χ² test statistic
    edges = binedges(fh)
    counts = bincounts(fh)
    n = length(data)
    expected_probs = CDF(edges[2:end]) .- CDF(edges[1:end-1])
    χsq = sum((counts .- n .* expected_probs) .^ 2 ./ (n .* expected_probs))
    return χsq
end

"""
    KolmogorovSmirnov(data, G)

Compute the Kolmogorov–Smirnov test statistic comparing empirical data to CDF `G`.

# Arguments
- `data`: Sample data.
- `CDF`:  Theoretical cumulative distribution function (CDF).

# Returns
The Kolmogorov–Smirnov statistic `√n * max|Fₙ(x) - G(x)|`.
"""
function KolmogorovSmirnov(data, CDF; m=1e-2)
    n = length(data)
    supp = collect(minimum(data):m:maximum(data))
    Fn = edf(data, supp)
    D = abs.(Fn .- CDF.(supp))
    return sqrt(n) * maximum(D)
end

"""
    AndersonDarling(data, G)

Compute the Anderson–Darling test statistic comparing `data` to CDF `G`.

# Arguments
- `data`: Sample data.
- `CDF`:  Theoretical cumulative distribution function (CDF).

# Returns
The Anderson–Darling statistic.
"""
function AndersonDarling(data, CDF)
    n = length(data)
    U = CDF.(sort(data))
    i = 1:n
    S = sum((2 .* i .- 1) .* (log.(U) .+ log.(1 .- reverse(U)))) / n
    return -(n + S)
end

"""
    CramervonMises(data, CDF)

Compute the Cramér–von Mises test statistic comparing `data` to CDF `G`.

# Arguments
- `data`: Sample data.
- `CDF`:  Theoretical cumulative distribution function (CDF).

# Returns
The Cramér–von Mises statistic.
"""
function CramervonMises(data, G)
    n = length(data)
    U = G.(sort(data))
    i = 1:n
    W = sum((U .- (2 .* i .- 1) ./ 2n) .^ 2)
    return W + 1 / (12 * n)
end

########################
### HELPER FUNCTIONS ###
"""
    edf

Construct the empirical distribution function (eDF/eCDF) from the data given some thresholds.
That is, computed F(x,t;n) = 1/n ∑ 1(x ≤ t), with 1(⋅) the indicator function.
"""
function edf(x::AbstractVector{<:Real}, t::AbstractVector{<:Real})
    n = length(x)
    F = [sum(sort(x) .<= τ) for τ in t] ./ n
    return F
end

end # module GOF
#/ End module
