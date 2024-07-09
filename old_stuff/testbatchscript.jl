# push!(LOAD_PATH, "/Users/millermacdonald/Desktop/Research_shit/NBI_Research")
# using dsnbIOmodule
include("dsnbIOmodule.jl")

using ArgParse


es_dsnb = range(0.5, 40, 100)

function params(p1, p2, p3)
    array = [p1, p2, p3]
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
#        "--case"
#            help = "sets branching ratios: 'B', 'C', or 'democratic'"
#            arg_type = String
        "--prog"
            help = "supernova progenitor: small, large, or black hole"
            arg_type = String
#        "--nubar"
#            help = "neutrino (true) or antineutrino (false)"
#            arg_type = Bool
        "--alpha1"
            help = "nu1 decay parameter"
            arg_type = Float64
        "--alpha2"
            help = "nu2 decay parameter"
            arg_type = Float64
    end
    parsed_args = parse_args(s)
    @time vals = params(parsed_args["alpha2"], parsed_args["alpha1"], parsed_args["prog"])
end

main()