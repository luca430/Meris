#= Module to plot Taylor's law for some processes and/or data =#
#/ Start module
module TaylorPlotter

using CairoMakie
using MakiePublication
using LaTeXStrings

using CSV, DataFrames, DataFramesMeta, JLD2

#################
### FUNCTIONS ###
function plot(
        filename;
        savefig = false,
        figname = "taylor.png"
    )

    @load filename taylor_d
    sc = Cycle([:color=>:markercolor, :strokecolor=>:color, :marker], covary=true)
    __theme = MakiePublication.theme_acs(; scattercycle=sc, ishollowmarkers=[true,true])
    set_theme!(__theme)

    colors = MakiePublication.COLORS[begin]

    width = .95 * 246
    height = 3*width / 4.67
    fig = Figure(; size=(1.5*width,height), figure_padding=(2,4,2,14))
    axtl = Axis(
        fig[1,1],
        xlabel=L"\textrm{sample mean}\;m", xlabelsize=12,
        ylabel=L"\textrm{sample variance}\;s^2", ylabelsize=12,
        # limits=(-15,0,-25,0)
    )

    #/ Plot
    xtl = -20:0.1:0.0
    classes = collect(keys(taylor_d))
    for class in classes[1:2]
        x, y = log.(taylor_d[class].means), log.(taylor_d[class].vars)
        scatter!(axtl, x, y, markersize=2, label=L"\text{%$class}")
    end

    for class in classes[1:2]
        b = taylor_d[class].fit.b
        ytl = b .* xtl
        lines!(axtl, xtl, ytl, linestyle=:dash, label=L"b = %$(round(b, digits=2))")
    end

    lines!(axtl, xtl, 2 .* xtl, color=:black)

    leg = Legend(fig[1, 2], axtl; tellheight=false)
    
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module TLPlotter
#/ End module
