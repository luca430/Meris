#= Simple module with formulae for pdfs and cdfs of Pareto-like distributions =#
#/ Start module
module ParetoLike

#################
### FUNCTIONS ###
function Paretocdf(α::Float64; xmin=0.0)
    F(x) = x < xmin ? 0.0 : 1 - (xmin / x)^α
    return F
end

function Burrcdf(c::Float64, α::Float64, λ::Float64)
    F(x) = 1 - (1 + (x / λ)^c)^(-α)
    return F
end

function Lomaxcdf(α::Float64, λ::Float64)
    return Burrcdf(1.0, α, λ)
end

end # module ParetoLike
#/ End module
