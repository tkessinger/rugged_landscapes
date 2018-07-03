#!/usr/bin/env julia

## sample_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic implementation of pop_sim

# tell Julia where the module is located
include("pop_sim.jl")
using Distributions, PopSim

K = 10000 # carrying capacity
fitnesses = [0, -0.01, 0.01, 0.1]
μ = 1e-4
# this mutation matrix allows 00 to mutate to 01, 10, or 11 at a reduced rate, etc.
# using a mutation matrix also means we can customize models very acutely, e.g., removing back mutations

numgens = 10000

# replace this with simple_forward_population(K,μ) to disallow back mutations
pop = simple_population(K, μ)

while pop.generation < numgens;
    evolve!(pop)
end

println("clones:")
[println(clone) for clone in pop.clones]
println("genotype frequencies")
println(get_frequencies(pop))