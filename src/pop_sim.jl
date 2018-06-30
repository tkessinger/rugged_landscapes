#!/usr/bin/env julia

## pop_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic simulation module on rugged fitness landscapes

using Distributions

## STILL TO DO
## population size and allele frequencies can be calculated more intelligently
## write some functions to generate a mutation matrix without the need to specify the whole thing
## implement some basic fitness landscape functions like HoC, RM, etc.
## consider making fitnesses and mutation matrix attributes of a "Landscape" type
## make the main functions class functions instead of constantly having to pass in the population

type Clone
    # new type for individual clones in population
    # each clone has a population size, genotype, and history
    # might be smart to make fitness an attribute as well, to avoid need to look up fitnesses constantly

    n::Int64
    genotype::Int64
    history::String

    # constructors
    Clone() = new(1,1,"1,0")
    Clone(n) = new(n,1,"1,0")
    Clone(n,genotype) = new(n,genotype,"1,0")
    Clone(n,genotype,history) = new(n,genotype,history)
end

function selection(population::Array{Any,1},K::Int64)
    # generates a poisson-distributed offspring number for each clone
    # note: presumably population wants to be population::Array{Clone, 1}--should this work?

    new_population = []
    mean_fitness = 0
    N = sum([clone.n for clone in population])

    # calculate the mean fitness
    # there might be a smarter way to do this other than iterating over all clones twice...
    for (cli, clone) in enumerate(population);
        mean_fitness += fitnesses[clone.genotype]*clone.n/N
    end

    # poisson offspring number distribution dependent on K and the fitness
    # the K/N term keeps the population size around K
    for (cli, clone) in enumerate(population);
        offspring_dist = Poisson(clone.n*(1+fitnesses[clone.genotype]-mean_fitness)*K/N)
        num_offspring = rand(offspring_dist)
        clone.n = num_offspring
        # prune any empty clones
        if clone.n != 0
            push!(new_population, clone)
        end
    end
    return new_population
end

function mutation(population::Array{Any,1},mut_matrix::Array{Float64,2}, mut_matrix_sum::Array{Float64,2},generation::Int64)
    # based on the sizes of each clone, sends some number of offspring to a new clone.
    # the i,j element of the mutation matrix is the mutation rate from genotype i to j (or should be)
    # mut_matrix_sum being passed in avoids the need to sum over the i axis repeatedly

    new_population = []
    for (cli, clone) in enumerate(population);
        # sample a number of individuals to make mutants
        # we will later decide what states they mutate into
        mut_dist = Binomial(clone.n,mut_matrix_sum[clone.genotype])
        num_muts = rand(mut_dist)
        if num_muts != 0
            clone.n -= num_muts
            # send the mutant individuals to different genotypes, weighted by their mutation rates
            mut_offspring = Multinomial(num_muts,mut_matrix[clone.genotype,:]/mut_matrix_sum[clone.genotype])
            new_clones_n = rand(mut_offspring)
            for (gi, ni) in enumerate(new_clones_n);
                # add a new clone if it's not empty
                if ni != 0
                    push!(new_population, Clone(ni,gi, string(clone.history,"::",gi,",",generation)))
                end
            end
        end
        # delete this clone if it's empty
        if clone.n != 0
            push!(new_population, clone)
        end
    end
    return new_population
end

function get_frequencies(population::Array{Any,1}, numloci::Int64)
    # spit out the genotype frequencies from the list of clones
    # there may be a more elegant way to do this
    freqs = zeros(numloci)
    N = sum([clone.n for clone in population])
    for (cli, clone) in enumerate(population);
        freqs[clone.genotype] += clone.n/N
    end
    return freqs
end

function cleanup(population::Array{Any,1})
    # idiotproofing function that removes empty clones
    for (cli, clone) in enumerate(population);
        if clone.n == 0
            splice!(population, cli)
        end
    end
    return population
end

K = 1000 # carrying capacity
numloci = 2 # number of loci
alph_length = 2 # length of the alphabet/number of states per locus
# note that we can also define a custom fitness landscape with no links whatsoever to an allele model
fitnesses = [0 -0.01 0.01 0.1]
mut_rate = 1e-2
# this mutation matrix allows 00 to mutate to 01, 10, or 11 at a reduced rate, etc.
# using a mutation matrix also means we can customize models very acutely, e.g., removing back mutations
mut_matrix = [0 mut_rate mut_rate mut_rate^2; mut_rate 0 mut_rate^2 mut_rate; mut_rate mut_rate^2 0 mut_rate; mut_rate^2 mut_rate mut_rate 0]
mut_matrix_sum = sum(mut_matrix,1)

numgens = 10000
generation = 0

# initialize the population
population = [Clone(K,1)]

while generation < numgens;
    population = selection(population, K)
    population = mutation_with_matrix(population, mut_matrix, mut_matrix_sum, generation)
    generation += 1
end

# get the state of the population by calling "population" or "get_frequencies(population, 4)"
