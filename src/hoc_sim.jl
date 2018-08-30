#!/usr/bin/env julia

## hoc.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulates a house-of-cards population until fixation.
## Parses the population's history and records it.

using Distributions, PopSim, ArgParse, JLD

# K = 100000 # carrying capacity
# numloci = 3
# numstates = 2
# μ = 1e-4
# steepness = 1e-3
#
# numtrials = 10
#
# outfile = "test"

# initialize the settings (the description is for the help screen)
function main(args)
    s = ArgParseSettings(description = "Example 1 for argparse.jl: minimal usage.")

    @add_arg_table s begin
        "--job_name"
            arg_type=AbstractString
            default = "test"
        "--K"
            arg_type=Float64
            default=1e5
        "--u"
            arg_type=Float64
            default=1e-5
        "--steepness"
            arg_type=Float64
            default=1e-2
        "--numloci"
            arg_type=Int64
            default=3
        "--numstates"
            arg_type=Int64
            default=2
        "--numtrials"
            arg_type=Int64
            default = 100
        "--file"
            arg_type=AbstractString
            default = "test"
    end

    parsed_args = parse_args(s) # the result is a Dict{String,Any}

    K = floor.(Int64,parsed_args["K"])

    job_name, μ, steepness, numloci, numstates, numtrials, outfile = parsed_args["job_name"], parsed_args["u"], parsed_args["steepness"], parsed_args["numloci"], parsed_args["numstates"], parsed_args["numtrials"], parsed_args["file"]

    println("running $numtrials simulations:")

    file = ismatch(r"\.jld", outfile) ? outfile : outfile*".jld"

    parsed_hists = []
    valley_percentage = []
    landscapes = []

    for trial in range(1,numtrials);
        println("initiating trial $trial")
        pop = hoc_population(K, numloci, numstates, μ, steepness)
        while !exists_path(numloci, numstates, pop.landscape.fitnesses);
            println("generating new landscape")
            pop = hoc_population(K, numloci, numstates, μ, steepness)
        end
        println(pop.landscape.fitnesses)
        evolve_until_fixation!(pop)
        parsed_hist, valleys, ridges = parse_history(pop)
        push!(parsed_hists, parsed_hist)
        push!(valley_percentage, valleys)
        push!(landscapes, pop.landscape.fitnesses)

        save("output/hoc_sims/$file", "hists", parsed_hists, "valleys", valley_percentage, "landscapes", landscapes)

    end
    save("output/hoc_sims/$file", "hists", parsed_hists, "valleys", valley_percentage, "landscapes", landscapes)
end

main(ARGS)
