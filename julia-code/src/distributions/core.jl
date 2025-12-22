#= Module to separate core and non-core components in a dataset =#
#/ Start module
module CoreFinder

#/ Packages
using Optim, StatsBase

function critical_point(x; slimits=(-5.0,8.0), n=1002, tol=1e-2)
    mm = moments(x; slimits=slimits, n=n)

    sgrid = mm.sgrid
    Λvals = mm.Λvals
    
    sgrid1 = mm.sgrid1
    Λprime = mm.Λprime
    
    sgrid2 = mm.sgrid2
    Λsecond = mm.Λsecond
    
    # Get critical point
    idx = findlast(Λsecond .> maximum(Λsecond) * tol)
    s_c = sgrid1[idx+1]
    x_c = Λprime[idx+1]

    return (x_c = x_c, s_c = s_c)
end

function moments(x; slimits=(-5.0,8.0), n=1002)
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
        sgrid1  = collect(sgrid[2:end-1]),
        Λprime   = -Λprime,
        sgrid2  = collect(sgrid[3:end-2]),
        Λsecond  = Λsecond
    )
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

end # module RFCSampler
#/ End module
