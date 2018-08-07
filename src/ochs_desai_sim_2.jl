#!/usr/bin/env julia

## ochs_desai_sim.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Simulate a series of populations as in Ochs and Desai (2015):
## competition between a simple selective sweep and crossing a fitness valley.
## Reproduces figure 3b.

using Distributions, PopSim, ArgParse, JLD

# arguments (with default values):
# K = Int64(100) carrying capacity
# slist = Float64[0.05, 0.0, 0.07]
# μlist = Float64[5.0*10.0^-6, 5.0*10.0^-5, 5.0*10.0^-5]
# 1 u -- rate from wild type to simple adaptation
# 2 i -- rate from wild type to valley
# 3 v -- rate from valley to complex adaptation
# genotype 2 is the simple adaptation and 4 is the complex one

function ochs_desai_ensemble(pop::Population, verbose=false)
    # quickly simulates an Ochs-Desai population
    # returns 0 if the sweep occurred and 1 if the valley was crossed
    # it might be smart to generalize this function
    freqs = get_frequencies(pop)
    while freqs[2] != 1 && freqs[4] != 1;
        evolve!(pop)
        freqs = get_frequencies(pop)
        if verbose && (pop.generation % 1000 == 0)
            println("generation $(pop.generation), N = $(sum([x.n for x in pop.clones])), frequencies = $freqs")
        end
    end
    if freqs[2] == 1
        # if the simple sweep has occurred
        return 0
    elseif freqs[4] == 1
        # if the valley was crossed
        return 1
    end
end

function compute_γ(slist::Array{Float64},μlist::Array{Float64})
    s_u, δ, s_v = slist
    u = μlist[1]
    if δ < 2*sqrt(u*s_v)
        #γ = 1.0/s_u*(0.5*log((s_u/u)^2) - (s_u/s_v)*log(s_u/u) + pi^2/6)
        γ = 1.0/s_u*(0.5*log(s_u/u)^2 - (s_u/s_v)*log(s_u/u) + pi^2/6)
    elseif δ > 2*sqrt(u*s_v)
        γ = 1.0/δ*(log(s_u/u) - (s_u/δ)*(1.0-(s_u/u)^(-δ/s_u)))
    end
    return γ
end

function main(args)

    # initialize the settings (the description is for the help screen)
    s = ArgParseSettings(description = "Example 1 for argparse.jl: minimal usage.")

    @add_arg_table s begin
        "--job_name"
            arg_type=AbstractString
            default = "test"
        "--NGamma"
            arg_type=Float64
            default=0.15
        "--s_u"
            arg_type=Float64
            default=0.0005
        "--s_v"
            arg_type=Float64
            default=0.05
        "--u"
            arg_type=Float64
            default=1e-6
        "--delta"
            arg_type=AbstractString
            default="high"

        "--numtrials"
            arg_type=Int64
            default = 100
        "--file"
            arg_type=AbstractString
            default = "test"
    end

    parsed_args = parse_args(s) # the result is a Dict{String,Any}

    u = parsed_args["u"]
    s_u = parsed_args["s_u"]
    s_v = parsed_args["s_v"]


    if parsed_args["delta"] == "high"
        δ = 10*sqrt(u*s_v)
    elseif parsed_args["delta"] == "low"
        δ = 0.0
    else
        println("invalid value of δ")
    end


    slist = [s_u, δ, s_v]
    μlist = [u, u, u]

    γ = compute_γ(slist,μlist)

    Γ = u^2*s_v*γ/s_u

    K = floor.(Int64,parsed_args["NGamma"]/Γ)

    println("K = $K, Γ = $Γ, γ = $γ, s = $slist, μ = $μlist")

    numtrials = parsed_args["numtrials"]
    outfile = parsed_args["file"]

    results = []

    params = OchsDesaiParams(K,slist,μlist)

    # 1 u -- rate from wild type to simple adaptation
    # 2 i -- rate from wild type to valley
    # 3 v -- rate from valley to complex adaptation

    #println("$parsed_args")

    #println("running $numtrials simulations:")

    for trial in range(1,numtrials);
        println("initiating trial number $trial")
        pop = ochs_desai_population(K, μlist, slist)
        push!(results,ochs_desai_ensemble(pop))
        file = ismatch(r"\.jld", outfile) ? outfile : outfile*".jld"
        save("output/ochs_desai_sims/$file", "params", params, "results", mean(results))
        # println(trial)
    end

    #println("proportion of runs where the valley was crossed:")
    #println(mean(results))

    file = ismatch(r"\.jld", outfile) ? outfile : outfile*".jld"
    save("output/ochs_desai_sims_2/$file", "params", params, "results", mean(results))

end

main(ARGS)
