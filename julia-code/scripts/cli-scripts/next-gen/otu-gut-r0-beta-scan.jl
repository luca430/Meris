#= One-dimensional scan for P(R0 > 1) over mean(beta) at selected biomass sigma_B. =#

using Meris

module R0Grid
include("otu-gut-r0-grid.jl")
end

const DEFAULT_OUTFILE_1D = joinpath(Meris.DATADIR, "next-gen", "otu-gut1-r0-beta-scan.jld2")
const DEFAULT_BIOMASS_MEAN_1D = "1000.0"
const DEFAULT_BIOMASS_SIGMAS_1D = "0.1,1.0,2.5"
const DEFAULT_N_GRID_1D = "50"

function print_help()
    println("""
    Usage:
      julia --project=julia-code --threads=10 julia-code/scripts/cli-scripts/next-gen/otu-gut-r0-beta-scan.jl [options]

    One-dimensional version of otu-gut-r0-grid.jl. It scans mean(beta)
    and uses three default sigma_B values at fixed mean biomass.

    Default overrides:
      --outfile=$(DEFAULT_OUTFILE_1D)
      --biomass-mean=$(DEFAULT_BIOMASS_MEAN_1D)
      --biomass-sigmas=$(DEFAULT_BIOMASS_SIGMAS_1D)
      --n-grid=$(DEFAULT_N_GRID_1D)

    All options accepted by otu-gut-r0-grid.jl are accepted here as well.
    User-supplied options override these one-dimensional defaults.
    """)
end

function main(args=ARGS)
    if "--help" in args || "-h" in args
        print_help()
        return nothing
    end

    default_args = [
        "--outfile=$(DEFAULT_OUTFILE_1D)",
        "--biomass-mean=$(DEFAULT_BIOMASS_MEAN_1D)",
        "--biomass-sigmas=$(DEFAULT_BIOMASS_SIGMAS_1D)",
        "--n-grid=$(DEFAULT_N_GRID_1D)",
    ]
    return R0Grid.main(vcat(default_args, args))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
