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
export fitness_ridge_population, fitness_ridge_population_mod, hoc_population

export OchsDesaiParams

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

    # μsum[i] is just the total mutation flux out of state i
    # this avoids the need to repeatedly compute the sum

    # constructor
    Landscape(fitnesses,μmatrix) = new(fitnesses, μmatrix, sum(μmatrix,2))
end

struct OchsDesaiParams
    K::Int64
    slist::Array{Float64,1}
    μlist::Array{Float64,1}
end

mutable struct Population
    # type for storing both the population itself (clone list) and associated parameters
    K::Int64
    landscape::Landscape
    clones::Array{Clone,1}
    generation::Int64

    # constructor
    Population(K,landscape,clones) = new(K,landscape,clones,0)
    Population(K,landscape) = new(K,landscape,Clone[Clone(K,1)],0)
end



function genotype_mapping(genotype::Int64, numloci::Int64, numstates::Int64)
    # maps a genotype number to an array containing the actual state of each locus
    genotype_list = Int64[]
    genotype -= 1
    for i in reverse(0:numloci-1);
        push!(genotype_list, div(genotype,numstates^i))
        genotype = genotype%(numstates^i)
    end
    return genotype_list
end



function greater_than_zero_condition(array1::Array{Int64}, array2::Array{Int64})
    # ensures that, anywhere array1 is nonzero, array2 is also nonzero
    conditionmet = true
    for (xi, x) in enumerate(array1);
        if array1[xi] != 0 && array2[xi] == 0
            conditionmet = false
        end
    end
    return conditionmet
end



function hoc_population(K::Int64, numloci::Int64, numstates::Int64, μ::Float64)
    fitnesses = hoc_fitnesses(numloci, numstates)
    μmatrix = hypercube_mutations(fitnesses, numloci, numstates,μ)
    landscape = Landscape(fitnesses,μmatrix)
    population = Population(K,landscape)
    return population
end



function hoc_fitnesses(numloci::Int64, numstates::Int64, steepness::Float64; α::Float64=0.0, first_constrained=true, final_highest=false)
    # generates random fitness values for the "house of cards" landscape model
    # α allows the "alpha-constrained" house of cards to be used,
    # where the first fitness is set to α
    # the last genotype fitness can be set to 1 with final_highest
    # setting first_constrained to false means the first fitness will be random

    # to do: add a function to allow every all-nonzero genotype to be a peak
    # i don't think this is necessary, but it might be a good idea

    numgenotypes = numstates^numloci

    fitnesses = Float64[]
    push!(fitnesses, first_constrained*α + (1-first_constrained)*rand())
    [push!(fitnesses, rand()) for x in 2:numgenotypes-1]
    push!(fitnesses, final_highest + (1-final_highest)*rand())

    fitnesses *= steepness

    return fitnesses
end


function rmf_fitnesses(numloci::Int64, numstates::Int64, θ::Float64, steepness::Float64)

    numgenotypes = numstates^numloci
    fitnesses = Float64[]
    norm = Normal(0,1)
    for x in 1:numgenotypes;
        dist = hamming_distance(genotype_mapping(x,numloci,numstates))
        rn = rand(norm)
        fit = θ*dist+rn
        push!(fitnesses,fit)
        println("$x, $dist, $θ, $rn,$fit")
    end
    println(fitnesses)
    fitnesses *= steepness
    return fitnesses

end


function hamming_distance(genotype1::Array{Int64}, genotype2::Array{Int64})
    # gives the Hamming distance between two genotypes
    return sum([genotype1[x] != genotype2[x] for x in 1:length(genotype1)])
end

function hamming_distance(genotype1::Array{Int64})
    # gives the Hamming distance between one genotype and the wild type (all zeros)
    return sum([genotype1[x] != 0 for x in 1:length(genotype1)])
end


function hypercube_mutations(fitnesses::Array{Float64}, numloci::Int64, numstates::Int64, μ::Float64;
    disallow_back::Bool=false, absorbing::Bool=true)
    # simple function for generating a "hypercube" mutation model
    # this can probably be typed/generalized: consider Hypercube_Landscape
    # currently this only works for 2 states per locus (numstates = 2)
    # we will generalize this later

    # disallow_back: forbid mutations back to state 0
    # absorbing: forbid mutations out of the max fitness state

    numgenotypes = numstates^numloci
    μmatrix = zeros(numgenotypes, numgenotypes)

    jumps = allowed_jumps(numloci, numstates)

    genotypelist = [genotype_mapping(n, numloci, numstates) for n in 1:numgenotypes]

    # figure out which loci are maximum fitness
    # if "absorbing" is on, we will prevent mutations out of the max fitness genotype

    if absorbing
        maxfit = maximum(fitnesses)
        maxlocs = find(x->x==maxfit,fitnesses)
    end

    # fill in the mutation matrix
    # if n loci need to be "flipped" during a particular jump,
    # then the rate needs to be set to μ^n
    # currently only 0 to non-0 flips are allowed

    for (fi, from_genotype) in enumerate(1:numgenotypes);
        # if the genotype is max fitness and absorbing is on, prevent out-mutations
        if absorbing && fi in maxlocs
            μmatrix[from_genotype,:] = 0
        # otherwise, proceed as normal
        else
            from_list = genotypelist[from_genotype]
            for (ti, to_genotype) in enumerate(1:numgenotypes);
                to_list = genotypelist[to_genotype]
                # if disallow_back: disallow jumps from, e.g., 01 to 10, but allow jumps from 01 to 02
                # essentially, jumps between nonzero states at each locus are allowed
                # mutation rate from a genotype to itself is of course zero
                if (greater_than_zero_condition(from_list, to_list) && from_list != to_list) || !disallow_back
                    boolsum = sum([from_list[x] != to_list[x] for x in 1:numloci])
                    if boolsum > 0
                        μmatrix[from_genotype, to_genotype] = μ^boolsum
                    end
                end
            end
        end
    end

    # note that if we wanted to allow back mutations, we could just transpose the mutation matrix
    # assuming, of course, that the rates are symmetric

    #landscape = Landscape(fitnesses,μmatrix)

    #clones = Clone[Clone(K,1)]

    #return Population(K, landscape, clones)
    return μmatrix
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

    return Population(K, landscape)
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

    return Population(K, landscape)

end

function simple_forward_population(K::Int64, μrate::Float64)
    # initializes a simple population with 2 loci and a valley
    # forward mutations only
    # state 4 is the final state: state 2 is a ridge crossing, 3 is a valley
    fitnesses = [0,0.01,-0.01,0.1]
    μmatrix = [0 μrate μrate μrate^2; 0 0 0 μrate; 0 0 0 μrate; 0 0 0 0]

    landscape = Landscape(fitnesses,μmatrix)

    return Population(K, landscape)
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
