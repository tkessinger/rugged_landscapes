#!/usr/bin/env julia

## ochs_desai_submit.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Submit script for ochs_desai_sims.jl;
## runs Ochs-Desai sims in batch using SLURM.

using JLD

pars = ["K", "s_u", "s_i", "s_v", "u_u", "u_i", "u_v"]
Klist = floor.(Int64,collect(logspace(2,6,100)))
sulist = [.006,.012,.025,.051]
parvals = [Klist, sulist, [0.0], [0.07], [5e-6], [5e-5], [5e-5]]
parsets = collect(Base.product(parvals...))
numtrials = 1000

nsets = length(parsets)

basename = "ochs_desai_sims"

for p in 1:nsets

    simstr = "src/ochs_desai_sim.jl"

    # create filname base from basename and run number
    numstr = lpad(string(p), length(string(nsets)), "0")
    filebase = basename * "_" * numstr

    # add parameter values as command line options
    simstr *= " "
#    simstr *= join(map( (x) -> "--" * x[1] * "=\"" * string(x[2]) * "\"", zip(pars, parsets[p])), " ")
# simstr *= " --numtrials=\"$numtrials\" "
# simstr *= " --file=\"" * filebase * ".jld\""

    simstr *= join(map( (x) -> "--" * x[1] * " " * string(x[2]), zip(pars, parsets[p])), " ")
    simstr *= " --numtrials $numtrials"
    simstr *= " --file " * filebase * ".jld"

    # cmdstr = "--job_name='ochs_desai_" * numstr * "' " *
    #     "--output=" * basename * ".out" *
    #     " --wrap='" * simstr * "'"

    #cmdstr = "--wrap=\'" * simstr * "\'"

#    print("sbatch src/ochs_desai_sim.jl $cmdstr")

#    run(`sbatch src/ochs_desai_sim.jl $cmdstr`)
    print("sbatch $(split(simstr))\n")
    run(`sbatch $(split(simstr))`)
end
