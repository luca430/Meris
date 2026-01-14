module Taylor

import Meris
using DataFrames, DataFramesMeta, StatsBase
using JLD2

function compute(
        df;
        nbins::Int = 27,
        save = false,
        filename = "tl.jld2"
    )

    #/ Get specific taylor df. Do this one class at a time to preserve memory in the computation
    classes = unique(df.class)
    taylor_d = Dict()
    for class in classes
        println(class)
        sub = df[df.class .== class, :]
        tldf = Meris.Taylor.compute(sub, :component_id, minoccupancy=1e-2, maxfrequency=1.0, occ=false)
    
        m = tldf[!,:meanfrequency]
        s = tldf[!,:varfrequency]
        x = log.(m)
        y = log.(s)

        #~ Compute rescaled moments using the occupancy        
        tldf = @chain tldf begin 
            @transform(:omeanfrequency = :meanfrequency .* :occupancy)
            @transform(:ovarfrequency = :varfrequency .+ (1 .- :occupancy) .* :meanfrequency.^2)
            @transform(:ovarfrequency = :ovarfrequency .* :occupancy)
        end
        
        #/ Fix a straight line using York's method
        #/ Calculate weights using the errors
        #~ Calculate how much of the total variation comes from presence-absence
        #  recall (σ′)² ← o⋅[σ²+μ²(1-o)], and so the ratio R = (σ′)² / (o⋅σ²) 
        o = tldf[!,:occupancy] .* (1 .- tldf[!,:occupancy])
        R = o .* tldf[!,:meanfrequency].^2 ./ tldf[!,:ovarfrequency]
        #~ filter those with ratio 0 [occupancy 0]
        sidxs = findall(x -> x > 0, R)
        x = x[sidxs]
        y = y[sidxs]
        #~ extract the errors on the mean and variance [see `taylor.jl`]
        #! note: use the δ-method to get the error on the log-transformed variables
        σx = m[sidxs] ./ m[sidxs].^2
        σy = s[sidxs] ./ s[sidxs].^2
        logcov = tldf[!,:errorcov][sidxs] ./ (m[sidxs] .* s[sidxs])
        logρ = logcov ./ sqrt.(σx .* σy)
        #~ specify the weights
        #! note: As for the line fitting only the relative weights are relevant, one could in
        #        principle scale the weights such that they are numerically more 'stable'. Yet,
        #        this may distort the error on the slope and intercept, as these are now
        #        'artificially' inflated by the weights. To bring them into a reasonable scale,
        #        we here specify the scale specifically, such that errors are reflecting the
        #        actual scatter of the means and variances and not the artificial weights.
        wx = 1.0 ./ σx
        wy = 1.0 .* sqrt.(1.0 .- R[sidxs]) ./ σy
        wscale = length(sidxs) / sum(wy)
        wx = wx .* wscale
        wy = wy .* wscale
        #~ Fit
        fit = Meris.StraightLine.weightedyorkfit(x, y, wx, wy, ρ=logρ)
        means = tldf[!, :omeanfrequency]
        vars = tldf[!, :ovarfrequency]
        taylor_d[class] = (means = means, vars = vars, fit = fit)
    end

    (save) && (@save filename taylor_d)

    return taylor_d
end

end # End module