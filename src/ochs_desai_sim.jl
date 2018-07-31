#!/usr/bin/env julia

## ochs_desai_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulate a series of populations as in Ochs and Desai (2015):
## competition between a simple selective sweep and crossing a fitness valley.

# tell Julia where the module is located
# include("pop_sim.jl")
using Distributions, PopSim, ArgParse, JLD

# arguments (with default values):
# K = Int64(100) carrying capacity
# slist = Float64[0.05, 0.0, 0.07]
# μlist = Float64[5.0*10.0^-6, 5.0*10.0^-5, 5.0*10.0^-5]
# 1 u -- rate from wild type to simple adaptation
# 2 i -- rate from wild type to valley
# 3 v -- rate from valley to complex adaptation
# genotype 2 is the simple adaptation and 4 is the complex one

function ochs_desai_ensemble(pop)
    # quickly simulates an Ochs-Desai population
    # returns 0 if the sweep occurred and 1 if the valley was crossed
    # it might be smart to generalize this function
    freqs = get_frequencies(pop)
    while freqs[2] != 1 && freqs[4] != 1;
        evolve!(pop)
        freqs = get_frequencies(pop)
    end
    if freqs[2] == 1
        # if the simple sweep has occurred
        return 0
    elseif freqs[4] == 1
        # if the valley was crossed
        return 1
    end
end


function main(args)

    # initialize the settings (the description is for the help screen)
    s = ArgParseSettings(description = "Example 1 for argparse.jl: minimal usage.")

    @add_arg_table s begin
        "--job_name"
            arg_type=AbstractString
            default = "test"
        "--K"   # carrying capacity
            arg_type=Float64
            default = 100.0
        "--s_u"   # list of selection coefficients
            arg_type=Float64
            default = 0.05
        "--s_i"
            arg_type=Float64
            default = 0.0
        "--s_v"
            arg_type=Float64
            default = 0.07
        "--u_u"   # list of mutation rates
            arg_type=Float64
            default = 5e-6
        "--u_i"   # list of mutation rates
            arg_type=Float64
            default = 5e-5
        "--u_v"   # list of mutation rates
            arg_type=Float64
            default = 5e-5
        "--numtrials"
            arg_type=Int64
            default = 100
        "--file"
            arg_type=AbstractString
            default = "test"
    end

    parsed_args = parse_args(s) # the result is a Dict{String,Any}

    K = floor.(Int64,parsed_args["K"])
    slist = [parsed_args["s_u"], parsed_args["s_i"], parsed_args["s_v"]]
    μlist = [parsed_args["u_u"], parsed_args["u_i"], parsed_args["u_v"]]

    numtrials = parsed_args["numtrials"]
    outfile = parsed_args["file"]

    results = []

    params = OchsDesaiParams(K,slist,μlist)

    # 1 u -- rate from wild type to simple adaptation
    # 2 i -- rate from wild type to valley
    # 3 v -- rate from valley to complex adaptation

    #println("$parsed_args")

    #println("running $numtrials simulations:")

    for trial in range(1,numtrials);
        pop = ochs_desai_population(K, μlist, slist)
        push!(results,ochs_desai_ensemble(pop))
        # println(trial)
    end

    #println("proportion of runs where the valley was crossed:")
    #println(mean(results))

    file = ismatch(r"\.jld", outfile) ? outfile : outfile*".jld"
    save("output/ochs_desai_sims/$file", "params", params, "results", mean(results))

end

main(ARGS)
