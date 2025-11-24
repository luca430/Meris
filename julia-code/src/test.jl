#/ Start module
module Test

using Distributions
using Random
using StatsBase

#################
### FUNCTIONS ###
import Meris.ParetoDistribution as MPareto

function sampleandfit(; α=0.5, β=2.0, τ=1e3, ε=1.0, n::Int=10^5)
    #~ Notes
    #  Zipf's scaling with exponent -1 implies β=2
    #  We typically should have β>α
    #  τ is the "inflection" point, or the typical scale at which another scaling "takes over"
    #  x > ε, so ε == xmin
    dPareto = MPareto.DoublePareto(α, β, τ, ε)
    rng = Random.Xoshiro(42)
    X = MPareto.rand(rng, dPareto, n)
    #~ Fit using MLE
    #  note: fitting assumes some fixed ε
    dPareto_mle = MPareto.fit(MPareto.DoublePareto, X; ε=ε)
    @info "params" MPareto.params(dPareto) MPareto.params(dPareto_mle)
    return nothing
end


########################
### HELPER FUNCTIONS ###
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

function doublepowerlaw(; α::Float64=.5, β::Float64=1., τ::Float64=1e3, ε::Float64=1.0)
    function F(x,α,β,τ,ε)
        numer = (ε/τ)^(1/α) + (ε/τ)^(1/β)
        denum = (x/τ)^(1/α) + (x/τ)^(1/β)
        return numer / denum
    end
    function f(x,α,β,τ,ε)
        εv = (ε/τ)^(1/α) + (ε/τ)^(1/β)
        numer = εv * (α*(x/τ)^(1/β) + β*(x/τ)^(1/α))
        denum = α*β*x*((x/τ)^(1/β) + β*(x/τ)^(1/α))^2
        return numer / denum
    end
    xmin = log10(ε)
    xmax = xmin + 8
    x = exp10.(range(xmin, xmax, 128))
    Fy = F.(x, Ref(α), Ref(β), Ref(τ), Ref(ε))
    fy = x.^(4α) .* f.(x, Ref(α), Ref(β), Ref(τ), Ref(ε))
    return (; x=x, Fy=Fy, fy=fy)
end

end # module Test

#/ End module
