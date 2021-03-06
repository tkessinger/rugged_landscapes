#!/usr/bin/env julia

## ridge_valley_submit.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Submit script for hoc_sim.jl;
## runs ridge-valley sims in batch using SLURM.

using JLD

pars = ["K", "s_u", "s_v", "u_u", "u_v", "ridgelength"]
Klist = floor.(Int64,collect(logspace(2,6,10)))
sulist = [0.01, 0.1]
svlist = [0.001, 0.01]
uulist = [1e-4]
uvlist = [1e-3]
#ridgelengthlist = [2]
parvals = [Klist, sulist, svlist, uulist, uvlist]#, ridgelengthlist]
# take the Cartesian product of all parameter combinations
parsets = collect(Base.product(parvals...))
numtrials = 100
pars = ["K", "s_u", "s_v", "u_u", "u_v", "ridgelength"]

nsets = length(parsets)

basename = "ridge_valley_sims"

for p in 1:nsets

    simstr = "src/ridge_valley_sim.jl"

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
