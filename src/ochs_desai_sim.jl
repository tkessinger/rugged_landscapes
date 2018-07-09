#!/usr/bin/env julia

## ochs_desai_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulate a series of populations as in Ochs and Desai (2015):
## competition between a simple selective sweep and crossing a fitness valley.

# tell Julia where the module is located
include("pop_sim.jl")
using Distributions, PopSim

K = 100 # carrying capacity

slist = Float64[0.05, 0.0, 0.07]
μlist = Float64[5.0*10.0^-6, 5.0*10.0^-5, 5.0*10.0^-5]
# 1 u -- rate from wild type to simple adaptation
# 2 i -- rate from wild type to valley
# 3 v -- rate from valley to complex adaptation
# genotype 2 is the simple adaptation and 4 is the complex one

numtrials = 100

results = []

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

println("running $numtrials simulations:")

for trial in range(1,numtrials);
    pop = ochs_desai_population(K, μlist, slist)
    push!(results,ochs_desai_ensemble(pop))
    println(trial)
end

println("proportion of runs where the valley was crossed:")
println(mean(results))
