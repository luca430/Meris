#= Module to investigate macro(-ecological) laws in RFC document =#
#/ Start module
module CoreFinder

#/ Packages
using Optim, Interpolations, StatsBase

function criticalities(x; slimits=(-5.0,8.0), n=1002, effort=100)
    boundary = s_boundary(x; slimits=slimits, n=n, effort=effort)
    s_c = boundary.s_c
    
    Λvals = boundary.Λvals
    sgrid = boundary.sgrid
    
    Λprime_s = boundary.Λprime_s
    Λprime = boundary.Λprime
    
    Λsecond_s = boundary.Λsecond_s
    Λsecond = boundary.Λsecond
    
    # interpolate Λ'(s_c) from the grid
    itp = LinearInterpolation(Λprime_s, Λprime)
    Λprime_at_sc = itp(s_c)
    
    # canonical expectation x(s) = -Λ'(s)
    x_c = -Λprime_at_sc
    n_c = exp(-x_c) * length(x)

    return (s_c = s_c, x_c = x_c, n_c = n_c)
end


#### HELPER FUNCTIONS ####
# -----------------------------------------------------------------------------
# 1. Stable log-mean-exp for z = vector. This trick avoids overflow.
# -----------------------------------------------------------------------------
logmeanexp(z::AbstractVector) = begin
    m = maximum(z)
    m + log(mean(exp.(z .- m)))
end

# -----------------------------------------------------------------------------
# 2. Λ(s, x) = log mean( exp( -s * x ) )
#    derivative: Λ'(s) = -E_s[x]
# -----------------------------------------------------------------------------
Λ(s, x) = logmeanexp(-s .* x)

# -----------------------------------------------------------------------------
# 3. Compute Λ, Λ', Λ'' on a uniform grid using central finite differences
# -----------------------------------------------------------------------------
function compute_derivs(x; slimits=(-5.0,5.0), n=1001)
    sgrid = range(slimits[1], slimits[2], length=n)
    Λvals = [Λ(s, x) for s in sgrid]

    Δs = step(sgrid)
    n1 = n - 2
    n2 = n - 4

    # first derivative: central differences, aligned with sgrid[2:end-1]
    Λprime = similar(Λvals, n1)
    for i in 2:n-1
        Λprime[i-1] = (Λvals[i+1] - Λvals[i-1]) / (2Δs)
    end

    # second derivative: central differences, aligned with sgrid[3:end-2]
    Λsecond = similar(Λvals, n2)
    for i in 3:n-2
        Λsecond[i-2] = (Λvals[i+1] - 2Λvals[i] + Λvals[i-1]) / (Δs^2)
    end

    return (
        sgrid    = sgrid,
        Λvals    = Λvals,
        s_cent1  = collect(sgrid[2:end-1]),
        Λprime   = Λprime,
        s_cent2  = collect(sgrid[3:end-2]),
        Λsecond  = Λsecond
    )
end

# -----------------------------------------------------------------------------
# 4. Find s_c considering a specific effort
# -----------------------------------------------------------------------------
function s_boundary(x; slimits=(-5.0,5.0), n=1001, effort=100)
    D = compute_derivs(x; slimits=slimits, n=n)
    # idx = argmax(D.Λsecond)
    idx = findlast(D.Λsecond .> maximum(D.Λsecond) / effort)
    s_c = D.s_cent2[idx]
    return (
        s_c = s_c,
        sgrid = D.sgrid,
        Λvals = D.Λvals,
        Λprime_s = D.s_cent1,
        Λprime = D.Λprime,
        Λsecond_s = D.s_cent2,
        Λsecond = D.Λsecond
    )
end


end # module RFCSampler
#/ End module
