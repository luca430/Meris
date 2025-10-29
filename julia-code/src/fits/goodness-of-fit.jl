#= Module to perform goodness-of-fit evaluations to check how well some data
   - is described by a particular distribution
   - is described by some distribution compared to another distribution
=#
#/ Start module
module GOF

#/ Packages
using Random
using StatsBase
using FHist

#################
### FUNCTIONS ###
"""
    estimatep(data, CDF, invCDF, test; nruns=10000, rng=Random.Xoshiro(42), params...)

Estimate p-value for a goodness-of-fit (GOF) test using Monte-Carlo simulations.

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
function edf(x::Array{Float64}, t::Array{Float64})
    n = length(x)
    F = [sum(sort(x) .<= τ) for τ in t] ./ n
    return F
end

end # module GOF
#/ End module
