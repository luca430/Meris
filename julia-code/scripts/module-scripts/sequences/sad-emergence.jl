#= SAD generator from heavy-tailed latent variables =#
#/ Start module
module SADGenerator

#/ Packages
using Distributions
using StatsBase
using Random
using FHist
using JLD2

import Meris.DATADIR as DATADIR
import Meris.MDistributions as MDist

#################
### FUNCTIONS ###
function generate(
    N::Int;
    K::Int=1,
    S::Int=10^4,
    γ::Float64=0.5,
    φ::Float64=1e8,
    ε::Float64=1.0,
    fitdistribution=false,
    verbose=true
)
    rng = Random.Xoshiro(N*42)
    Pθ = MDist.TemperedPareto(γ, 1/φ, ε)

    (verbose) && (println("Sampling..."))
    n = zeros(Int, S, K)
    p = zeros(Float64, S, K)
    intensity = 0.
    
    for k in 1:K
        θ = MDist.rand(rng, Pθ, S)
        intensity += sum(θ)
        p[:,k] = sort(θ) ./ sum(θ)
        Mult = Distributions.Multinomial(N, p[:,k])
        counts = rand(rng, Mult)
        n[:,k] .+= counts
    end
    n = reduce(vcat, n)
    nbar = n[n.>0] / (N*K)
    
    #/ Compute histogram
    (verbose) && (println("Fitting..."))
    logn = log10.(nbar)
    bmin, bmax = extrema(logn)
    fhbinedges = range(bmin, bmax, 31)
    fh = FHist.Hist1D(logn; binedges=fhbinedges, overflow=true) |> normalize
    p = fh.bincounts[fh.bincounts .> 0]
    n = bincenters(fh)[fh.bincounts .> 0]

    !(fitdistribution) && (return (; n=n, p=p))
    
    #/ Fit a (tempered) Pareto distribution
    nmin, nmax = extrema(logn)
    εs = exp10.(range(nmin, nmax/100, 63))
    # Pfit = MDist.fit(MDist.TemperedPareto, nbar; εs=εs)
    Pfit = MDist.TemperedPareto(γ, intensity/φ, exp10(-5.75))
    paretox = exp10.(range(log10(Pfit.ε+1e-9), 1+nmax, 127))
    Z = sum(count(x -> x > Pfit.ε, nbar)) / length(nbar)
    paretoy = log(10) .* paretox .* MDist.pdf.(Pfit, paretox) * Z
    
    #~ Fit a normal to the "bulk"
    # nbulk = nbar[nbar .<= Pfit.ε]
    # Norm = Distributions.fit_mle(Distributions.Normal, nbar)
    # normx = exp10.(range(log10(minimum(nbar)), log10(Pfit.ε), 63))
    # normy = log(10) .* normx .* Distributions.pdf.(Norm, normx)
    
    return (; n=n, p=p, paretox=paretox, paretoy=paretoy)
end

function save_sad(result, filename; DIR=DATADIR*"sad/synthetic/", fitdistribution=true)
	  #~ Store the synthetic SAD and its params
    mkpath(DIR)
    if fitdistribution
        jldsave(
            DIR*filename;
            n=result.n, p=result.p,
            paretox=result.paretox, paretoy=result.paretoy
        )
        return nothing
    end
    jldsave(DIR*filename; n=result.n, p=result.p)
    return nothing
end

end # module SADGenerator
#/ End module
