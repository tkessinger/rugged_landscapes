#!/usr/bin/env julia

## ridge_valley_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulate competition between a ridge-traversing and valley-crossing
## path to a particular complex adaptation.

using Distributions, PopSim, JLD, ArgParse

function ridge_valley_ensemble(pop::Population)
    # simulates a ridge/valley crossing competition
    # returns 0 if the ridge was crossed, 1 if the valley was crossed, and 2 if indeterminate
    freqs = get_frequencies(pop)
    ridgelength = length(freqs)-3
    while freqs[ridgelength+1] != 1 && freqs[ridgelength+3] != 1;
        evolve!(pop)
        freqs = get_frequencies(pop)
    end
    if freqs[ridgelength+1] == 1
        # if the ridge was traversed
        return 0
    elseif freqs[ridgelength+3] == 1
        # if the valley was crossed
        return 1
    else
        println("error: ridge_valley_ensemble failed to fix either adaptation")
        return 2
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
            default = 0.1
        "--s_v"
            arg_type=Float64
            default = 0.01
        "--u_u"   # list of mutation rates
            arg_type=Float64
            default = 1e-6
        "--u_v"   # list of mutation rates
            arg_type=Float64
            default = 1e-5
        "--ridgelength"
            arg_type=Int64
            default = 2
        "--numtrials"
            arg_type=Int64
            default = 10
        "--file"
            arg_type=AbstractString
            default = "test"
    end

    parsed_args = parse_args(s) # the result is a Dict{String,Any}

    K = floor.(Int64,parsed_args["K"])
    slist = [parsed_args["s_u"], -parsed_args["s_v"]]
    μlist = [parsed_args["u_u"], parsed_args["u_v"]]
    ridgelength = parsed_args["ridgelength"]

    numtrials = parsed_args["numtrials"]
    outfile = parsed_args["file"]

    file = ismatch(r"\.jld", outfile) ? outfile : outfile*".jld"

    # K = 100000 # carrying capacity
    # slist = Float64[0.1,-0.01]
    # μlist = Float64[10.0^-3,10.0^-2]
    # ridgelength = 2
    # 1 ridge mutation rate, also height of ridge (s)
    # 2 valley mutation rate, also depth of valley (δ--should be negative)
    # intermediate steps across the ridge will be set to s/n, with n the ridgelength
    # genotype ridgelength+1 will be the complex adaptation
    # genotype ridgelength+2 will be the valley state

    #numtrials = 100

    results = []


    println("running $numtrials simulations:")

    for trial in range(1,numtrials);
        pop = fitness_ridge_population(K, μlist, slist, 2)
        push!(results,ridge_valley_ensemble(pop))
        println("initializing $trial")
        save("output/ridge_valley_sims/$file", "results", mean(results), "numruns", trial)
    end
    save("output/ridge_valley_sims/$file", "results", mean(results), "numruns", numtrials)

    println("proportion of runs where the valley was crossed:")
    println(mean(results))
end

main(ARGS)
