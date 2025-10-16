module GOF

using Statistics, Random, FHist

"""
    make_EDF(data, x)

Compute the empirical distribution function (EDF) of `data` evaluated at points `x`.

# Arguments
- `data`: Sample data.
- `x`: Points at which to evaluate the EDF.

# Returns
A vector of EDF values `Fₙ(x) = (1/n) * sum(data ≤ x)`.
"""
function make_EDF(data, x)
    n = length(data)
    sorted_data = sort(data)
    return [sum(sorted_data .<= xi) / n for xi in x]
end


"""
    chi2(data, G; nbins=30, nparams=0)

Perform a χ² (chi-squared) goodness-of-fit test comparing `data` to distribution `G`.

# Arguments
- `data`: Sample data.
- `G`: Cumulative distribution function (CDF) of the theoretical distribution.
- `nbins`: Number of histogram bins (default = 30).
- `nparams`: Number of estimated parameters (optional, currently unused).

# Returns
The chi-squared test statistic.
"""
function chi2(data, G; nbins=30, nparams=0)
    bmin, bmax = round(minimum(data)), round(maximum(data))
    Δb = (bmax - bmin) / nbins
    fh = FHist.Hist1D(data, binedges=bmin:Δb:bmax)
    edges = binedges(fh)
    counts = bincounts(fh)

    n = length(data)
    expected_probs = G(edges[2:end]) .- G(edges[1:end-1])
    X2 = sum((counts .- n .* expected_probs) .^ 2 ./ (n .* expected_probs))
    return X2
end


"""
    KS(data, G; supp=-10:1e-3:10)

Compute the Kolmogorov–Smirnov test statistic comparing empirical data to CDF `G`.

# Arguments
- `data`: Sample data.
- `G`: Theoretical cumulative distribution function (CDF).
- `supp`: Support grid for evaluating the EDF (default = -10:1e-3:10).

# Returns
The Kolmogorov–Smirnov statistic `√n * max|Fₙ(x) - G(x)|`.
"""
function KS(data, G; supp=-10:1e-3:10)
    n = length(data)
    Fn = make_EDF(data, supp)
    D = abs.(Fn .- G.(supp))
    return sqrt(n) * maximum(D)
end


"""
    AD(data, G)

Compute the Anderson–Darling test statistic comparing `data` to CDF `G`.

# Arguments
- `data`: Sample data.
- `G`: Theoretical CDF.

# Returns
The Anderson–Darling statistic.
"""
function AD(data, G)
    n = length(data)
    U = G.(sort(data))
    i = 1:n
    A = sum((2 .* i .- 1) .* (log.(U) .+ log.(1 - reverse(U))))
    return -A / n - n
end


"""
    CvM(data, G)

Compute the Cramér–von Mises test statistic comparing `data` to CDF `G`.

# Arguments
- `data`: Sample data.
- `G`: Theoretical CDF.

# Returns
The Cramér–von Mises statistic.
"""
function CvM(data, G)
    n = length(data)
    U = G.(sort(data))
    i = 1:n
    W = sum((U .- (2 .* i .- 1) ./ n) .^ 2)
    return W + 1 / (12n)
end


"""
    GOF(data, G, G_inv, test; nruns=10000, rng=Random.Xoshiro(42), params...)

Monte Carlo (MC) estimation of the p-value for a goodness-of-fit (GOF) test.

# Arguments
- `data`: Sample data.
- `G`: Theoretical cumulative distribution function (CDF).
- `G_inv`: Inverse CDF (quantile function) used for sampling under the null hypothesis.
- `test`: A function implementing the test statistic (e.g., `KS`, `AD`, `CvM`, `chi2`).
- `nruns`: Number of Monte Carlo runs (default = 10,000).
- `rng`: Random number generator (default = `Random.Xoshiro(42)`).
- `params...`: Additional keyword arguments passed to the test function.

# Returns
Estimated p-value `p = mean(T_sim ≥ T)` where `T` is the observed statistic.
"""
function MC_gof(data, G, G_inv, test; nruns=10000, rng=Random.Xoshiro(42), params...)
    n = length(data)
    T = test(data, G; params...)
    samples = [G_inv.(rand(rng, Float64, n)) for _ in 1:nruns]
    T_sim = [test(s, G; params...) for s in samples]
    p = mean(T_sim .>= T)
    return p
end

end # module Utils
