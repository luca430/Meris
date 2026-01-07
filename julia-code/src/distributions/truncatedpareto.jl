struct TruncatedPareto{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    ε::T
    εmax::T
end
