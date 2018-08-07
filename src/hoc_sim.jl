#!/usr/bin/env julia

## hoc.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulates a house-of-cards population until fixation.
## Parses the population's history and records it.

using Distributions, PopSim, JLD

K = 100000 # carrying capacity
numloci = 4
numstates = 2
μ = 1e-4
steepness = 1e-5

numtrials = 100

outfile = "test"

println("running $numtrials simulations:")

file = ismatch(r"\.jld", outfile) ? outfile : outfile*".jld"

parsed_hists = []
valley_percentage = []
landscapes = []

for trial in range(1,numtrials);
    println("initiating trial $trial")
    pop = hoc_population(K, numloci, numstates, μ, steepness)
    println(pop.landscape.fitnesses)
    evolve_until_fixation!(pop)
    parsed_hist, valleys, ridges = parse_history(pop)
    push!(parsed_hists, parsed_hist)
    push!(valley_percentage, valleys)
    push!(landscapes, pop.landscape.fitnesses)

    save("output/$file", "hists", parsed_hists, "valleys", valley_percentage, "landscapes", landscapes)

end
