#!/usr/bin/env julia

## pop_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic simulation module on rugged fitness landscapes

## STILL TO DO
## population size and allele frequencies can be calculated more intelligently
## write some functions to generate a mutation matrix without the need to specify the whole thing
## implement some basic fitness landscape functions like HoC, RM, etc.
## consider making fitnesses and mutation matrix attributes of a "Landscape" type

## bug: back mutation rate seems to spike very late in simulations
## either that, or else history is not working correctly

module PopSim
export evolve!, get_frequencies, simple_forward_population, simple_population

using Distributions

struct Clone
    # new type for individual clones in population
    # each clone has a population size, genotype, and history
    # might be smart to make fitness an attribute as well, to avoid need to look up fitnesses constantly
    # struct is NOT mutable, be wary. we might change this
    # history is an array of 1x2 arrays: first element is the genotype, second is the mutation time

    n::Int64
    genotype::Int64
    history::Array{Array{Int64,1}}

    # constructors
    Clone() = new(1,1,[[1,0]])
    Clone(n) = new(n,1,[[1,0]])
    Clone(n,genotype) = new(n,genotype,[[1,0]])
    Clone(n,genotype,history) = new(n,genotype,history)
end

struct Landscape
    # type for storing fitness function and mutation matrix
    # i will probably write a function to allow explicit genotypes to be stored, too

    fitnesses::Array{Float64,1}
    μmatrix::Array{Float64,2}
    μsum::Array{Float64,2}

    Landscape(fitnesses,μmatrix) = new(fitnesses, μmatrix, sum(μmatrix,2))
end

mutable struct Population
    # type for storing both the population itself (clone list) and associated parameters
    K::Int64
    landscape::Landscape
    clones::Array{Clone,1}
    generation::Int64

    # constructors
    Population(K,landscape,clones) = new(K,landscape,clones,0)
end

function simple_forward_population(K::Int64, μrate::Float64)
    # initializes a simple population with 2 loci and a valley
    # forward mutations only
    fitnesses = [0,0.01,-0.01,0.1]
    μmatrix = [0 μrate μrate μrate^2; 0 0 0 μrate; 0 0 0 μrate; 0 0 0 0]

    landscape = Landscape(fitnesses,μmatrix)
    clones = Clone[Clone(K,1)]
    return Population(K, landscape, clones)
end

function simple_population(K::Int64, μrate::Float64)
    # initializes a simple population with 2 loci and a valley
    # forward mutations only
    fitnesses = [0,0.01,-0.01,0.1]
    μmatrix = [0 μrate μrate μrate^2; μrate 0 μrate^2 μrate; μrate μrate^2 0 μrate; μrate^2 μrate μrate 0]

    # create landscape and clone list
    landscape = Landscape(fitnesses,μmatrix)
    clones = Clone[Clone(K,1)]
    return Population(K, landscape, clones)
end


function change_clone_size(clone::Clone, newsize::Int64)
    # copies an existing clone but modifies the size
    new_clone = Clone(newsize, clone.genotype, clone.history)
    return new_clone
end

function selection!(population::Population)
    # generates a poisson-distributed offspring number for each clone
    # note: presumably population wants to be population::Array{Clone, 1}--should this work?

    new_clones = Clone[]
    mean_fitness = 0
    N = sum([clone.n for clone in population.clones])

    # calculate the mean fitness
    # there might be a smarter way to do this other than iterating over all clones twice...
    for (cli, clone) in enumerate(population.clones);
        mean_fitness += population.landscape.fitnesses[clone.genotype]*clone.n/N
    end

    # poisson offspring number distribution dependent on K and the fitness
    # the K/N term keeps the population size around K
    for (cli, clone) in enumerate(population.clones);
        offspring_dist = Poisson(clone.n*(1+population.landscape.fitnesses[clone.genotype]-mean_fitness)*population.K/N)
        num_offspring = rand(offspring_dist)
        #clone.n = num_offspring

        tmp_clone = change_clone_size(clone, num_offspring)
        # prune any empty clones
        if tmp_clone.n != 0
            push!(new_clones, tmp_clone)
        end
    end
    # update clone list
    population.clones = new_clones
end

function mutation!(population::Population)
    # based on the sizes of each clone, sends some number of offspring to a new clone.
    # the i,j element of the mutation matrix is the mutation rate from genotype i to j
    # mut_matrix_sum being passed in avoids the need to sum over the i axis repeatedly

    new_clones = Clone[]
    for (cli, clone) in enumerate(population.clones);
        # sample a number of individuals to make mutants
        # we will later decide what states they mutate into
        # only call this block if the total mutation rate out of a genotype is > 0
        parent_clone_size = clone.n
        if population.landscape.μsum[clone.genotype] > 0
            mut_dist = Binomial(clone.n,population.landscape.μsum[clone.genotype])
            num_muts = rand(mut_dist)
            if num_muts > 0
                parent_clone_size -= num_muts
                # send the mutant individuals to different genotypes, weighted by their mutation rates
                mut_offspring = Multinomial(num_muts,population.landscape.μmatrix[clone.genotype,:]/population.landscape.μsum[clone.genotype])
                new_clones_n = rand(mut_offspring)
                for (gi, ni) in enumerate(new_clones_n);
                    # add a new clone if it's not empty
                    if ni != 0
                        newhist = copy(clone.history)
                        push!(newhist, [gi, population.generation])
                        push!(new_clones, Clone(ni, gi, newhist))
                    end
                end
            end
        end
        # delete this clone if it's empty
        if parent_clone_size != 0
            push!(new_clones, change_clone_size(clone, parent_clone_size))
        end
    end
    # update clone list
    population.clones = new_clones
end

function get_frequencies(population::Population)
    # spit out the genotype frequencies from the list of clones
    # there may be a more elegant way to do this
    freqs = zeros(length(population.landscape.fitnesses))
    N = sum([clone.n for clone in population.clones])
    for (cli, clone) in enumerate(population.clones);
        freqs[clone.genotype] += clone.n/N
    end
    return freqs
end

function cleanup!(population::Population)
    # idiotproofing function that removes empty clones
    for (cli, clone) in enumerate(population.clones);
        if clone.n == 0
            splice!(population.clones, cli)
        end
    end
end

function evolve!(population::Population)
    selection!(population)
    mutation!(population)
    population.generation += 1
end

end
