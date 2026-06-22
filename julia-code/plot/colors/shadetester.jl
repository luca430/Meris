#/ Start module
module Shades

using StatsBase
using CairoMakie, Colors

function shades(base::Colorant, n)
    hsl_base = convert(HSL, base)
    h, s = hsl_base.h, hsl_base.s
    _colors = [HSL(h, s, l) for l in range(0.25, 0.85, length=n)]
    return _colors
end

function mean_rgb(cols)
    rgbs = RGB.(cols)
    RGB(mean(c.r for c in rgbs),
        mean(c.g for c in rgbs),
        mean(c.b for c in rgbs))
end

function check(;
    bases = [
        colorant"#1f77b4",  # blue
        colorant"#ff7f0e",  # orange
        colorant"#9467bd",  # purple
        colorant"#2ca02c",  # green
        colorant"#d62728"   # red
    ],
    nshades = [7, 10, 5, 4, 5]
    )

    fig = Figure(; size=(246, 246))

    for (i, base) in enumerate(bases)
        ax = Axis(fig[div(i-1,3)+1, mod(i-1,3)+1])
        cols = shades(base, nshades[i])

        for (j, c) in enumerate(cols)
            scatter!(ax,
                     randn(40) .+ j,
                     randn(40) .+ j,
                     color = c,
                     markersize = 8
                     )
        end
    end

    fig
end

end # module Shades
#/ End module
