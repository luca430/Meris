#= Module that contains methods other than OLS for fitting a straight line =#
#/ Start module
module StraightLine

#/ Packages
using StatsBase
using Random

#################
### FUNCTIONS ###
"""
Fit straight line using the methods from York
Assumes
"""
function weightedyorkfit(X, Y, XWEIGHTS, YWEIGHTS;
    ρ=nothing,
    tol=1e-6,
    maxiterations::Int = 100
)
    #~ Initial guess of the parameters
    #  note: while one can use OLS, or any other regression, a simple estimation works as well
    b = 2.0
    
    #~ Instantiate and allocate
    (isnothing(ρ)) && (ρ = zeros(length(X)))
    α = sqrt.(XWEIGHTS .* YWEIGHTS)
    W = similar(X)
    B = similar(X)
    XBAR = 0.
    YBAR = 0.
    converged = false

    iterations = 0
    eps = Inf
    while eps > tol && iterations < maxiterations
        W = (XWEIGHTS .* YWEIGHTS) ./ (XWEIGHTS .+ b^2 .* YWEIGHTS .- (2*b).*ρ.*α)
        ZW = sum(W)
        XBAR = sum(W.*X) ./ ZW
        YBAR = sum(W.*Y) ./ ZW
        U = X .- XBAR
        V = Y .- YBAR
        B = W .* (U ./ YWEIGHTS .+ b .* V ./ XWEIGHTS - ρ.*(b.*U .+ V) ./ α)
        WB = W .* B
        bprime = sum(WB .* V) / sum(WB .* U)
        #~ Compute error as |b - b'|
        eps = abs(b .- bprime)
        b = bprime
        iterations += 1
    end
    converged = iterations < maxiterations
    #/ Compute a from b, X and Y
    ZW = sum(W)
    a = YBAR - b*XBAR
    #/ Compute standard errors
    #~ Compute least-squares-adjusted points
    x = XBAR .+ B
    #~ Compute σa and σb
    xbar = sum(W.*x) / ZW
    u = x .- xbar
    σb = sqrt(1 / sum(W .* u.^2))
    σa = sqrt(1 / ZW + (xbar*σb)^2)

    #/ Return
    return (; a=a, b=b, σa=σa, σb=σb)
end

end # module StraightLine
#/ End module
