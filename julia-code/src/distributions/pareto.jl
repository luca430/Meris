#= Simple module with formulae for pdfs and cdfs of Pareto-like distributions =#
#/ Start module
module ParetoLike

#################
### FUNCTIONS ###
# CDFs
function Paretocdf(α::Float64; xmin=0.0)
    F(x) = x < xmin ? 0.0 : 1 - (xmin / x)^α
    return F
end

function Burrcdf(c::Float64, α::Float64, λ::Float64)
    F(x) = 1 - (1 + (x / λ)^c)^(-α)
    return F
end

# PDFs
function Paretopdf(x::Array{T}, α::Float64; xmin=1.0) where T<:Real
    return Paretopdf.(x, Ref(α); xmin=xmin)
end

function Paretopdf(x::Float64, α::Float64; xmin=1.0)
    (x < xmin) && (return 0.0)
    return α*xmin^α * x^(-α-1)
end

function Burrpdf(x::Array{T}, c::Float64, α::Float64, λ::Float64) where T<:Real
	  return Burrpdf.(x, Ref(c), Ref(α), Ref(λ))
end

function Burrpdf(x::T, c::Float64, α::Float64, λ::Float64) where T<:Real
    (x <= 0) && (return 0.0)
	  return (c*α/λ) * (x/λ)^(c-1) * (1 + (x/λ)^c)^(-α-1)
end

function Lomaxcdf(α::Float64, λ::Float64)
    return Burrcdf(1.0, α, λ)
end

function Lomaxpdf(x::T, α::Float64, λ::Float64) where T<:Real
    return Burrpdf.(x, 1.0, Ref(α), Ref(λ))
end

end # module ParetoLike
#/ End module
