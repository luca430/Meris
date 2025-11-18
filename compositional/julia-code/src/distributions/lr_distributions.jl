module LRDistr

using SpecialFunctions, Statistics

function lr_gamma(z, α)
    return α * sqrt(trigamma(α)) .* z .+ α * digamma(α) .- exp.(z .* sqrt(trigamma(α)) .+ digamma(α)) .+ 0.5 * log(trigamma(α)) .- loggamma(α)
end

function lr_lognormal(z, σ)
    return -z .^ 2 ./ 2 .- log(sqrt(σ^2 * 2 * π))
end

function lr_betap(z, a, b)
    s = sqrt(trigamma(a) + trigamma(b))
    m = digamma(a) - digamma(b)
    B = gamma(a) * gamma(b) / gamma(a + b)
    return log(s / B) .+ a .* (s .* z .+ m) .- (a + b) .* log.(1 .+ exp.(s .* z .+ m))
end

end # module LRDistr
