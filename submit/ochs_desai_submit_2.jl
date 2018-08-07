#!/usr/bin/env julia

## ochs_desai_submit.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Submit script for ochs_desai_sims.jl;
## runs Ochs-Desai sims in batch using SLURM.

using JLD

pars = ["NGamma", "s_u", "s_v", "u", "delta"]
#NGammalist = collect(linspace(.01,.15,8))
#sulist = 5.0*collect(logspace(-4,-2,10))
NGammalist = collect(linspace(.01,.15,3))
sulist = 5.0*collect(logspace(-4,-2,3))
svlist = [0.05,0.01]
ulist = [1e-5,1e-6]
deltalist = ["low", "high"]
parvals = [NGammalist, sulist, svlist, ulist, deltalist]
# take the Cartesian product of all parameter combinations
parsets = collect(Base.product(parvals...))
numtrials = 1000

nsets = length(parsets)

basename = "ochs_desai_sims_2"

for p in 1:nsets

    simstr = "src/ochs_desai_sim_2.jl"

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
