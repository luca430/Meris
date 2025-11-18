#= Module to fit a continuous double power law and estimate the point where the behaviour changes =#
#/ Start module
module DoublePowerLaw

#/ Packages
using StatsBase, LsqFit, LinearAlgebra

using Meris

#############################
###### MAIN FUNCTIONS #######

# --- Fit function: takes ranks, freqs and returns fit object and derived results
function fit_broken_zipf(ranks, freqs; use_binning=false, nbins=40)
    # Optionally bin
    if use_binning
        xr, yr = log_bin(ranks, freqs, nbins=nbins)
        xdata = collect(xr)
        ydata = collect(yr)
    else
        xdata = collect(ranks)
        ydata = collect(freqs)
    end

    # Fit in log-space: target = log(y)
    ylog = log.(ydata)

    # initial guesses
    loga0 = log(maximum(ydata))
    α10 = 0.8      # first-slope guess (shallow)
    α20 = 1.6      # second-slope guess (steeper)
    rb0 = median(xdata)  # guess breakpoint ~ middle rank
    p0 = [loga0, α10, α20, log(rb0)]

    # perform nonlinear least squares
    fit = curve_fit(log_broken_pl, xdata, ylog, p0)

    # covariance estimate and standard errors
    J = fit.jacobian           # jacobian at solution
    res = ylog .- log_broken_pl(xdata, fit.param)
    n = length(ylog)
    m = length(fit.param)
    σ² = sum(res.^2) / (n - m)          # residual variance estimate
    cov = inv(J' * J) * σ²
    stderr = sqrt.(diag(cov))

    # derive readable parameter values
    loga, α1, α2, logr = fit.param
    se_loga, se_a1, se_a2, se_logr = stderr
    rb = exp(logr)
    se_rb = rb * se_logr   # delta-method: var(rb) ≈ (rb^2)*(var(logr))

    # frequency at break
    A = exp(loga)
    freq_rb = A * rb^(-α1)
    # estimate se for freq_rb via delta method approximately (neglecting covariances)
    # We'll compute approximate relative error using var(logA) and var(logr) and var(alpha1)
    var_loga = cov[1,1]
    var_a1 = cov[2,2]
    var_logr = cov[4,4]
    # log(freq_rb) = logA - α1*log(rb)
    var_logfreq_rb = var_loga + (log(rb))^2 * var_a1 + (α1^2) * var_logr
    se_freq_rb = exp(log(freq_rb)) * sqrt(var_logfreq_rb)

    return (
        fit = fit,
        param = fit.param,
        stderr = stderr,
        A = A, α1 = α1, α2 = α2, rb = rb,
        se = (se_loga, se_a1, se_a2, se_logr),
        se_rb = se_rb,
        freq_rb = freq_rb,
        se_freq_rb = se_freq_rb,
        cov = cov,
        xdata = xdata,
        ydata = ydata,
        residuals = res
    )
end

#############################
##### HELPER FUNCTIONS ######

# --- Log-binning to reduce scatter
# x and y are rank and frequency
function log_bin(x, y; nbins=40)
    xmin, xmax = minimum(x), maximum(x)
    edges = 10 .^ range(log10(xmin), log10(xmax), length=nbins+1)
    xb = Float64[]
    yb = Float64[]
    for i in 1:length(edges)-1
        idx = (x .>= edges[i]) .& (x .< edges[i+1])
        if any(idx)
            push!(xb, mean(x[idx]))
            push!(yb, mean(y[idx]))
        end
    end
    return xb, yb
end

# --- Broken power-law model in log-space
# params p = [loga, alpha1, alpha2, logrb]
# x is rank
function log_broken_pl(x, p)
    loga, α1, α2, logr = p
    rb = exp(logr)
    A = exp(loga)
    # create y in linear space then return log(y)
    y = similar(x, Float64)
    for (i, r) in enumerate(x)
        if r <= rb
            y[i] = A * r^(-α1)
        else
            y[i] = A * rb^(α2 - α1) * r^(-α2)
        end
    end
    return log.(y)
end

#############################

end # module MLEFit
#/ End module
