#= Simple module with formulae for pdfs and cdfs of Pareto-like distributions =#
#/ Start module
module ParetoLike

#################
### FUNCTIONS ###

#######################
### LOG LIKELIHOODS ###
"""
    loggeneralizedPareto(arg, μ, σ, ξ)

Log density function of the generalized Pareto distribution, with an expansion with ξ near zero.
"""
function loggeneralizedPareto(x::T, μ::T, σ, ξ) where T<:Real
    z = (x - μ) / σ
    expn = if abs(ξ) < 1e-5
        #~ Expansion for ξ near zero
        -z * (ξ + 1) + (z^2) * ξ * (ξ + 1) / 2 -
        (z^3) * (ξ^2) * (ξ + 1) / 3 +
        (z^4) * (ξ^3) * (ξ + 1) / 4
    else
        (-(1 + ξ) / ξ) * log(max(0, 1 + z * ξ))
    end
    return expn - log(σ)
end

############
### CDFs ###
function generalizedParetocdf(σ::Float64, ξ::Float64; xmin=0.0)
    function F(x)
        (x < xmin) && (return 0.0)
        z = (x - xmin) / σ
        if iszero(ξ)
            return 1 - exp(-z)
        end
        return 1 - (1 + ξ*z)^(-1/ξ)
    end
    return F
end

function Paretocdf(α::Float64; xmin=0.0)
    F(x) = x < xmin ? 0.0 : 1 - (xmin / x)^α
    return F
end

function Burrcdf(c::Float64, α::Float64, λ::Float64)
    F(x) = 1 - (1 + (x / λ)^c)^(-α)
    return F
end

############
### PDFS ###
#~ Generalized Pareto distribution
function generalizedParetopdf(x::Array{T}, σ::Float64, ξ::Float64; xmin=1.0) where T<:Real
	  return generalizedParetopdf.(x, Ref(σ), Ref(ξ); xmin=xmin)
end

function generalizedParetopdf(x::Float64, σ::Float64, ξ::Float64; xmin=1.0)
    (x < xmin) && (return 0.0)
    z = (x - xmin) / σ
    return (1 + ξ*z) ^ (-1 - 1 / ξ) / σ
end

#~ Pareto distribution
function Paretopdf(x::Array{T}, α::Float64; xmin=1.0) where T<:Real
    return Paretopdf.(x, Ref(α); xmin=xmin)
end

function Paretopdf(x::Float64, α::Float64; xmin=1.0)
    (x < xmin) && (return 0.0)
    return α*xmin^α * x^(-α-1)
end

#~ Burr distribution
function Burrpdf(x::Array{T}, c::Float64, α::Float64, λ::Float64) where T<:Real
	  return Burrpdf.(x, Ref(c), Ref(α), Ref(λ))
end

function Burrpdf(x::T, c::Float64, α::Float64, λ::Float64) where T<:Real
    (x <= 0) && (return 0.0)
	  return (c*α/λ) * (x/λ)^(c-1) * (1 + (x/λ)^c)^(-α-1)
end

#~ Lomax distribution
Lomaxcdf(α::Float64, λ::Float64) = Burrcdf(1.0, α, λ)

function Lomaxpdf(x::T, α::Float64, λ::Float64) where T<:Real
    return Burrpdf.(x, 1.0, Ref(α), Ref(λ))
end

end # module ParetoLike
#/ End module
