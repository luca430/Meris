module LRDistr

using SpecialFunctions, Statistics

include("./DataTools.jl")
using .DataTools

function lr_gamma(z, α)
    return α*sqrt(trigamma(α)) .* z .+ α*digamma(α) .- exp.(z .* sqrt(trigamma(α)) .+ digamma(α)) .+ 0.5*log(trigamma(α)) .- loggamma(α)
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

### Useful functions ###
function gamma_params(df; occ=0.99)
    freqs = DataTools.get_frequencies(df; occ = occ, rescale=false)
    μv = mean(freqs, dims=1)
    σv = std(freqs, dims=1)
    βv = μv .^ 2 ./ σv .^ 2
    return mean(βv), std(βv)
end

function betap_params(df; occ=0.99)
    freqs = DataTools.get_frequencies(df; occ = occ, rescale=false)
    μv = mean(freqs, dims=1)
    σv = std(freqs, dims=1)
    θ = mean(μv) # THIS IS AN ANSATZ
    βv = μv .^ 2 ./ σv .^ 2 .+ μv .* θ ./ σv .^ 2 .+ 2
    αv = μv .^ 3 ./ σv .^ 2 ./ θ .+ μv .^ 2 ./ σv .^ 2 .+ μv ./ θ
    return mean(αv), std(αv), mean(βv), std(βv)
end

end # module Utils
