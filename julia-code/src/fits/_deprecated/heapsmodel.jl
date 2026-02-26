#= Module that contains methods to fit three regimes Heaps model =#
#/ Start module
module HeapsModel

#/ Packages
using LsqFit
using Statistics

#################
### FUNCTIONS ###

# ---------- Models (log-log) ----------

# a < 1 : 3-regime model, p = [log10c1, Δ, log10m1, log10m2], with c2 = c1*10^Δ
function model_three_regime_log(logx, p; a::Float64)
    lc1, Δ, lm1, lm2 = p
    c1 = 10.0^lc1
    c2 = c1 * 10.0^Δ
    m1 = 10.0^lm1
    m2 = 10.0^lm2

    x = 10.0 .^ logx
    y = x .* (1 .+ (x ./ c1).^m1).^((a - 1)/m1) .*
            (1 .+ (x ./ c2).^m2).^(-a/m2)

    return log10.(y)
end

# a == 1 : reduced model, p = [log10c2, log10m2]
function model_a1_log(logx, p)
    lc2, lm2 = p
    c2 = 10.0^lc2
    m2 = 10.0^lm2

    x = 10.0 .^ logx
    y = x .* (1 .+ (x ./ c2).^m2).^(-1/m2)
    return log10.(y)
end

# ---------- Fit wrapper ----------

"""
    fit_regimes(xdata, ydata; a, a_tol=1e-12, bounds=:default, p0=nothing)

Fits the appropriate regime model for a given `a` with 0 < a <= 1.

- If `a` is (within `a_tol`) equal to 1, fits the reduced 2-regime model.
- Otherwise fits the full 3-regime model.

Returns `(fit, params_named, model_symbol)` where:
- `model_symbol` is `:a1` or `:three_regime`
- `params_named` contains parameters in linear space.
"""
function fit_regimes(xdata, ydata; a::Float64, a_tol=1e-12, p0=nothing,
                     m_bounds=(0.2, 50.0), Δ_bounds=(0.5, 20.0), c_bounds=(-30.0, 30.0))

    @assert 0.0 < a <= 1.0 "Require 0 < a <= 1"
    @assert all(xdata .> 0) "xdata must be > 0 for log-fit"
    @assert all(ydata .> 0) "ydata must be > 0 for log-fit (drop zeros or add pseudocount)"

    logx = log10.(xdata)
    logy = log10.(ydata)

    # Helper for bounds in log-parameter space
    (m_lo, m_hi) = m_bounds
    (Δ_lo, Δ_hi) = Δ_bounds
    (lc_lo, lc_hi) = c_bounds

    if abs(a - 1.0) <= a_tol
        # ---- a = 1 reduced model ----
        if p0 === nothing
            c2_0 = quantile(xdata, 0.8)
            m2_0 = 2.0
            p0 = [log10(c2_0), log10(m2_0)]
        end

        lower = [lc_lo, log10(m_lo)]
        upper = [lc_hi, log10(m_hi)]

        fit = curve_fit(model_a1_log, logx, logy, p0; lower=lower, upper=upper)

        lc2, lm2 = fit.param
        c2 = 10.0^lc2
        m2 = 10.0^lm2
        return fit, (a=1.0, c2=c2, m2=m2), :a1
    else
        # ---- a < 1 full 3-regime model ----
        if p0 === nothing
            c1_0 = quantile(xdata, 0.2)
            c2_0 = quantile(xdata, 0.8)
            Δ0   = max(log10(c2_0 / c1_0), 1.0)
            m1_0 = 2.0
            m2_0 = 2.0
            p0 = [log10(c1_0), Δ0, log10(m1_0), log10(m2_0)]
        end

        lower = [lc_lo, Δ_lo, log10(m_lo), log10(m_lo)]
        upper = [lc_hi, Δ_hi, log10(m_hi), log10(m_hi)]

        fit = curve_fit((lx,p)->model_three_regime_log(lx, p; a=a),
                        logx, logy, p0; lower=lower, upper=upper)

        lc1, Δ, lm1, lm2 = fit.param
        c1 = 10.0^lc1
        c2 = c1 * 10.0^Δ
        m1 = 10.0^lm1
        m2 = 10.0^lm2
        plateau = c1^(1-a) * c2^a

        return fit, (a=a, c1=c1, c2=c2, m1=m1, m2=m2, plateau=plateau), :three_regime
    end
end

# ---------- Optional: evaluate fitted curve in linear space ----------

function predict_regimes(x, params; a_tol=1e-12)
    if haskey(params, :c1)
        a  = params.a
        c1 = params.c1
        c2 = params.c2
        m1 = params.m1
        m2 = params.m2
        return x .* (1 .+ (x ./ c1).^m1).^((a - 1)/m1) .*
                    (1 .+ (x ./ c2).^m2).^(-a/m2)
    else
        # a == 1 reduced model
        c2 = params.c2
        m2 = params.m2
        return x .* (1 .+ (x ./ c2).^m2).^(-1/m2)
    end
end


end # module StraightLine
#/ End module
