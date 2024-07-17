
# push!(LOAD_PATH, "/Users/millermacdonald/Desktop/Research_shit/NBI_Research")
# using dsnbIOmodule
include("dsnb2numodule.jl")

using ArgParse

function params(p1, p2, p3)
    array = [p1, p2, p3]
end

es_dsnb = range(0.5, 40, 100)

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--channel"
            help = "sets branching ratios: 'B', 'C', or 'democratic'"
        "--prog"
            help = "supernova progenitor: small, large, or black hole"
            arg_type = String
#        "--nubar"
#            help = "neutrino (true) or antineutrino (false)"
#            arg_type = Bool
        "--alpha"
            help = "nui decay parameter"
            arg_type = Float64
    end
    parsed_args = parse_args(s)

    @time println([params(parsed_args["channel"], parsed_args["prog"], parsed_args["alpha"]),
    [DSNB_vdecay_2ν_νe_1pc_tot.(es_dsnb, parsed_args["alpha"], parsed_args["channel"], bar, "NO", parsed_args["prog"], SNRnorm) for bar in [true, false]]])
end

main()