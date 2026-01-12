#= SAD generator from heavy-tailed latent variables =#
#/ Start module
module SADGenerator

#/ Packages
using Distributions
using StatsBase
using Random

import Meris.MDistributions as MDist

#################
### FUNCTIONS ###
function generate(
    N::Int; K::Int=100, S::Int=10^4, γ::Float64=0.42, φ::Float64=1e8, ε::Float64=1.
)
    rng = Random.Xoshiro(N*42)
    Pθ = MDist.TemperedPareto(γ, 1/φ, ε)

    x = zeros(Float64, S)

    for k in 1:K
        θ = MDist.rand(rng, Pθ, S)
        p = θ ./ sum(θ)
        Mult = Distributions.Multinomial(N, p)
        counts = rand(rng, Mult)
        x .+= counts ./ N
    end
    return x / K
end

end # module SADGenerator
#/ End module
