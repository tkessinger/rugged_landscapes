#!/usr/bin/env julia

## pop_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic simulation module on rugged fitness landscapes

## STILL TO DO
## population size and allele frequencies can be calculated more intelligently
## write some functions to generate a mutation matrix without the need to specify the whole thing
## implement some basic fitness landscape functions like HoC, RM, etc.
## better handling of hypercube fitness landscapes
## map genotype number to string of loci values

module PopSim
export evolve!
export get_frequencies
export simple_forward_population, simple_population, ochs_desai_population
export fitness_ridge_population, fitness_ridge_population_mod

using Distributions, Combinatorics

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

struct Hypercube_Landscape
    # INCOMPLETE
    # type for storing fitness function and mutation matrix
    # this one includes an explicit genotype-"state" map
    # need to either subtype this into Landscape or vice versa


    fitnesses::Array{Float64,1}
    μmatrix::Array{Float64,2}
    μsum::Array{Float64,2}

    numloci::Int64
    numstates::Int64

    Hypercube_Landscape(fitnesses,μmatrix) = new(fitnesses, μmatrix, sum(μmatrix,2))
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

function single_jumps(numloci::Int64, numstates::Int64)
    # produces a list of all allowable single jump lengths,
    # for a given number of loci and alphabet length
    # e.g.: 2 loci, 2 states - single jumps are from 1 to 2 or 3,
    # or from 2 or 3 to 4:
    # allowable single jump sizes are 1 or 2.
    # this is currently broken insofar as we do not "remember"
    # which loci have already been "flipped".
    # see comment under hypercube_population.

    singlejumps = Int64[]
    for n in 0:numstates-2;
        for x in 1:numloci;
            push!(singlejumps,x*numstates^n)
        end
    end
    return singlejumps
end

function allowed_jumps(numloci::Int64, numstates::Int64)
    # produces a list of all allowable jumps of any size,
    # for each permitted size,
    # for a given number of loci and alphabet length
    singlejumps = single_jumps(numloci,numstates)
    jumps = [[sum(x) for x in combinations(singlejumps,n)] for n in 1:numloci]
    return jumps
end

function genotype_mapping(numloci::Int64, numstates::Int64)

end

function hypercube_population(K::Int64, numloci::Int64, numstates::Int64, μ::Float64)
    # simple function for generating a population with a "hypercube" fitness landscape
    # fitnesses are currently left empty, but the mutation rate matrix is populated
    # this can probably be typed/generalized: consider Hypercube_Landscape
    # currently this only works for 2 states per locus (numstates = 2)
    # we will generalize this later
    # this still needs to be debugged: the correct implementation is likely to be uglier
    # see note below

    numgenotypes = numstates^numloci
    fitnesses = zeros(numgenotypes)
    μmatrix = zeros(numgenotypes, numgenotypes)

    jumps = allowed_jumps(numloci, numstates)

    # fill in the mutation matrix
    # this is the tricky part.
    # we have a list of permitted jump sizes:
    # if n loci need to be "flipped" during a particular jump, then the rate needs to
    # be set to μ^n.

    # first loop ends at numgenotypes-1 because the end state is absorbing
    for (fi, from_genotype) in enumerate(1:numgenotypes-1);
        # second loop begins at from_genotype+1 because a state never mutates to itself
        for (ti, to_genotype) in enumerate(from_genotype+1:numgenotypes);
            for n in 1:numloci;
                if (to_genotype - from_genotype) in jumps[n]
                    μmatrix[from_genotype,to_genotype] = μ^n
                end
            end
        end
    end

    # note that if we wanted to allow back mutations, we could just transpose the mutation matrix
    # assuming, of course, that the rates are symmetric

    landscape = Landscape(fitnesses,μmatrix)

    clones = Clone[Clone(K,1)]

    return Population(K, landscape, clones)
end

function fitness_ridge_population(K::Int64, μlist::Array{Float64,1}, slist::Array{Float64,1}, ridgelength::Int64)
    # generates a "ridge competition" landscape
    # by convention μ and s should be as follows:
    # 1 ridge mutation rate, also height of ridge (s)
    # 2 valley mutation rate, also depth of valley (δ)
    # intermediate steps across the ridge will be set to s/ridgelength, times the step number
    # note that ridgelength does not include the wild type--the "real" length is ridgelength+1
    # ridgelength is more like the number of edges, not vertices
    # this is a "forward" model only: back mutations are disallowed

    ridge_fitnesses = Float64[x*slist[1]/ridgelength for x in 1:ridgelength]
    valley_fitnesses = Float64[slist[2]]
    fitnesses = vcat(Float64[0.0], ridge_fitnesses, valley_fitnesses)
    μmatrix = zeros(ridgelength+2,ridgelength+2)

    # this loop fills in the mutation matrix and allows for multiple mutations to occur
    # the rate from x to x+n should be μ^n
    # there should be ridgelength-n+1 such jumps
    for n in 1:ridgelength;
        [μmatrix[x,x+n] = μlist[1]^n for x in 1:ridgelength-n+1]
    end

    # set the valley crossing mutation rates
    μmatrix[1,ridgelength+2] = μlist[2]
    μmatrix[ridgelength+2,ridgelength+1] = μlist[2]
    # allow for a double mutation
    μmatrix[1,ridgelength+1] += μlist[2]^2
    landscape = Landscape(fitnesses,μmatrix)
    clones = Clone[Clone(K,1)]

    return Population(K, landscape, clones)

end

function fitness_ridge_population_mod(K::Int64, μlist::Array{Float64,1}, slist::Array{Float64,1}, ridgelength::Int64)
    # generates a "ridge competition" landscape
    # by convention μ and s should be as follows:
    # 1 ridge mutation rate, also height of ridge (s)
    # 2 valley mutation rate, also depth of valley (δ)
    # intermediate steps across the ridge will be set to s/ridgelength, times the step number
    # note that ridgelength does not include the wild type--the "real" length is ridgelength+1
    # ridgelength is more like the number of edges, not vertices
    # this is a "forward" model only: back mutations are disallowed
    # this may be a smarter implementation than fitness_ridge_population
    ridge_fitnesses = Float64[x*slist[1]/ridgelength for x in 1:ridgelength]
    valley_fitnesses = Float64[slist[2]]
    fitnesses = vcat(Float64[0.0], ridge_fitnesses, valley_fitnesses, slist[1])
    μmatrix = zeros(ridgelength+3,ridgelength+3)

    # this loop fills in the mutation matrix and allows for multiple mutations to occur
    # the rate from x to x+n should be μ^n
    # there should be ridgelength-n+1 such jumps
    for n in 1:ridgelength;
        [μmatrix[x,x+n] = μlist[1]^n for x in 1:ridgelength-n+1]
    end

    # set the valley crossing mutation rates
    μmatrix[1,ridgelength+2] = μlist[2]
    μmatrix[ridgelength+2,ridgelength+3] = μlist[2]
    # allow for a double mutation
    μmatrix[1,ridgelength+3] = μlist[2]^2
    landscape = Landscape(fitnesses,μmatrix)
    clones = Clone[Clone(K,1)]

    return Population(K, landscape, clones)
end


function ochs_desai_population(K::Int64, μlist::Array{Float64,1}, slist::Array{Float64,1})
    # generates a "competition" landscape qua Ochs and Desai (2015)
    # by convention μ and s should be as follows:
    # 1 u -- rate from wild type to simple adaptation
    # 2 i -- rate from wild type to valley
    # 3 v -- rate from valley to complex adaptation
    # genotype 2 is the simple adaptation and 4 is the complex one
    fitnesses = vcat(Float64[0.0], slist)
    μmatrix = Float64[0 μlist[1] μlist[2] μlist[2]*μlist[3]; 0 0 0 0; 0 0 0 μlist[3]; 0 0 0 0]

    landscape = Landscape(fitnesses,μmatrix)
    clones = Clone[Clone(K,1)]

    return Population(K, landscape, clones)

end

function simple_forward_population(K::Int64, μrate::Float64)
    # initializes a simple population with 2 loci and a valley
    # forward mutations only
    # state 4 is the final state: state 2 is a ridge crossing, 3 is a valley
    fitnesses = [0,0.01,-0.01,0.1]
    μmatrix = [0 μrate μrate μrate^2; 0 0 0 μrate; 0 0 0 μrate; 0 0 0 0]

    landscape = Landscape(fitnesses,μmatrix)
    clones = Clone[Clone(K,1)]
    return Population(K, landscape, clones)
end

function simple_population(K::Int64, μrate::Float64)
    # initializes a simple population with 2 loci and a valley
    # back mutations allowed
    # state 4 is the final state: state 2 is a ridge crossing, 3 is a valley
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

# final end statement to close the module
end
