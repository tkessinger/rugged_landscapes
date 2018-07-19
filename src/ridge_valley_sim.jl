#!/usr/bin/env julia

## ridge_valley_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulate competition between a ridge-traversing and valley-crossing
## path to a particular complex adaptation.

# tell Julia where the module is located
include("pop_sim.jl")
using Distributions, PopSim

K = 100000 # carrying capacity

# note that we are using fitness_ridge_population_mod in pop_sim
# this has slightly different behavior from vanilla fitness_ridge_population
slist = Float64[0.1,-0.01]
μlist = Float64[10.0^-3,10.0^-2]
ridgelength = 2
# 1 ridge mutation rate, also height of ridge (s)
# 2 valley mutation rate, also depth of valley (δ--should be negative)
# intermediate steps across the ridge will be set to s/n, with n the ridgelength
# genotype ridgelength+1 will be the complex adaptation
# genotype ridgelength+2 will be the valley state

numtrials = 100

results = []

function ridge_valley_ensemble(pop)
    # simulates a ridge/valley crossing competition
    # returns 0 if the ridge was crossed, 1 if the valley was crossed, and 2 if indeterminate
    freqs = get_frequencies(pop)
    ridgelength = length(freqs)-3
    while freqs[ridgelength+1] != 1 && freqs[ridgelength+3] != 1;
        evolve!(pop)
        freqs = get_frequencies(pop)
    end
    if freqs[ridgelength+1] == 1
        return 0
    elseif freqs[ridgelength+3] == 1
        # if the valley was crossed
        return 1
    else
        println("error: ridge_valley_ensemble failed to fix either adaptation")
        return 2
    end
end

println("running $numtrials simulations:")

for trial in range(1,numtrials);
    pop = fitness_ridge_population(K, μlist, slist, 2)
    push!(results,ridge_valley_ensemble(pop))
    println(trial)
end

println("proportion of runs where the valley was crossed:")
println(mean(results))
