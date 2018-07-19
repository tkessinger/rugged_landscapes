#!/usr/bin/env julia

## ochs_desai_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulate a series of populations as in Ochs and Desai (2015):
## competition between a simple selective sweep and crossing a fitness valley.

# tell Julia where the module is located
include("pop_sim.jl")
using Distributions, PopSim, ArgParse

function main(args)

    # initialize the settings (the description is for the help screen)
    s = ArgParseSettings(description = "Example 1 for argparse.jl: minimal usage.")

    @add_arg_table s begin
        "--K"                 # a positional argument
            arg_type=Int64
            default=100
        "--s"
            nargs=3
            arg_type=Float64
            default = [0.07, 0.0, 0.05]
        "--u"
            nargs=3
            arg_type=Float64
            default = [5e-6, 5e-5, 5e-5]
    end

    parsed_args = parse_args(s) # the result is a Dict{String,Any}
    println("Parsed args:")
    for (key,val) in parsed_args
        println("  $key  =>  $(repr(val))")
    end
end

main(ARGS)
