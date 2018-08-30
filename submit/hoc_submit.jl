#!/usr/bin/env julia

## hoc_submit.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Submit script for hoc_sim.jl;
## runs "house of cards sims in batch using SLURM.

using JLD

pars = ["K", "u", "steepness"]
#NGammalist = collect(linspace(.01,.15,8))
#sulist = 5.0*collect(logspace(-4,-2,10))
Klist = collect(logspace(3,6,7))
ulist = collect(logspace(-5,-2,5))
steepnesslist = collect(logspace(-2,-1,5))
numlocilist = [3,4]
parvals = [Klist, ulist, steepnesslist, numlocilist]
# take the Cartesian product of all parameter combinations
parsets = collect(Base.product(parvals...))
numtrials = 100

nsets = length(parsets)

basename = "hoc_sims"

for p in 1:nsets

    simstr = "src/hoc_sim.jl"

    # create filname base from basename and run number
    numstr = lpad(string(p), length(string(nsets)), "0")
    filebase = basename * "_" * numstr

    # add parameter values as command line options
    simstr *= " "
    simstr *= join(map( (x) -> "--" * x[1] * " " * string(x[2]), zip(pars, parsets[p])), " ")
    simstr *= " --numtrials $numtrials"
    simstr *= " --file " * filebase * ".jld"

    # run the scripts via sbatch
    # the "split" wizardry is needed due to the fine distinction between Cmd and String types
    print("sbatch $(split(simstr))\n")
    run(`sbatch $(split(simstr))`)
end
