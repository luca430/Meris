#/ Start module
module Test

using Distributions
using Random
using StatsBase

function f(; α::Float64 = 1.5, n::Int = 1000)
	  rng = Random.Xoshiro(42)
    p = Distributions.Pareto(α)
    X = rand(rng, p, n)

    #~ Zipf's
    r = (tiedrank(X, rev=true) .- 1) ./ n
    #~ Survival
    C = _ecdf(X)
    S = 1.0 .- C.F

    return (; r=r, x=X), (; S=S, x=C.x)
end

function _ecdf(xs::Array{T}; sorted=false) where T<:Real
    (!sorted) && (xs = sort(xs))
    n = length(xs)
    F = similar(xs, Float64)
    k = 1
    for i in eachindex(xs)
        #~ Move k until xs[k] > edges[k]
        while k ≤ n && xs[k] ≤ xs[i]
            k += 1
        end
        F[i] = (k-1) / n
    end
    return (; F=F, x=xs)
end

function doublepowerlaw(; α::Float64 = 1.5, τ::Float64 = 1e4, ε::Float64 = 1.0)
    function f(x,α,τ,ε)
        numer = (x + α*τ*(x/τ)^α) * (ε + τ*(ε/τ)^α)
        denum = x*(x + τ*(x/τ)^α)^2
        return numer / denum
    end
    xmin = log10(ε)
    xmax = xmin + 8
    x = exp10.(range(xmin, xmax, 128))
    y = x.^2 .* f.(x, Ref(α), Ref(τ), Ref(ε))
    return (; x=x, y=y)
end

end # module Test

#/ End module
