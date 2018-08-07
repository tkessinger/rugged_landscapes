#!/usr/bin/env julia

## PopSimBasics.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Basic implementation of PopSim types and functions.
## Mostly for testing purposes,

module PopSimBasics

export init_ridge_pop, init_simple_pop

using PopSim
using Distributions, Combinatorics

function init_ridge_pop()
    K = Int64(1e5)
    μlist = Float64[1e-5,1e-4]
    slist = Float64[0.01,-0.001]
    ridgelength = 4
    return fitness_ridge_population(K,μlist,slist,ridgelength)
end

function init_simple_pop()
    K = Int64(1e5)
    μrate = 1e-5
    return simple_forward_population(K,μrate)
end

function init_ochs_desai_pop()
    # generates a "competition" landscape qua Ochs and Desai (2015)
    # by convention μ and s should be as follows:
    # 1 u -- rate from wild type to simple adaptation
    # 2 i -- rate from wild type to valley
    # 3 v -- rate from valley to complex adaptation
    # genotype 2 is the simple adaptation and 4 is the complex one
    K = 1e5
    μlist = [1e-5,1e-4,1e-4]
    slist = [0.05, 0, 0.07]
    return ochs_desai_pop(K, μlist, slist)
end


end
