#!/usr/bin/env julia

## plt_ochs_desai.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Open outputs of ochs_desai_sim_2.jl and plot results.
## Duplicates plot 3b from Ochs and Desai (2015).

using PyPlot, JLD, Glob, PopSim, LaTeXStrings

# glob all the output files
files = glob("output/ochs_desai_sims_2/ochs_desai_sims_2_*.jld")

# keys for this dict will be NΓ, su, δ, sv, and μ
NΓ_results = Dict{Tuple{Float64, Float64, Float64, Float64, Float64}, Float64}()

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

function compute_Γ(slist::Array{Float64},μlist::Array{Float64})
    s_u, δ, s_v = slist
    u = μlist[1]
    γ = compute_γ(slist,μlist)
    Γ = u^2*s_v*γ/s_u
    return Γ
end

#γ = compute_γ(slist,μlist)

#Γ = u^2*s_v*γ/s_u

# open files and fill dict
for (fi, file) in enumerate(files)
    f = load(file)
    results = f["results"]
    params = f["params"]
    Γ = compute_Γ(params.slist,params.μlist)
    s_u, δ, s_v = params.slist
    μ = params.μlist[1]
    NΓ_results[params.K*Γ,s_u,δ,s_v,μ] = results
end

# pull out the unique key values
NΓ_vals = sort(unique(x[1] for x in keys(NΓ_results)))
su_vals = sort(unique(x[2] for x in keys(NΓ_results)))
δ_vals = sort(unique(x[3] for x in keys(NΓ_results)))
sv_vals = sort(unique(x[4] for x in keys(NΓ_results)))
μ_vals = sort(unique(x[5] for x in keys(NΓ_results)))

parsets = collect(Base.product([su_vals, δ_vals, sv_vals, μ_vals]...))

fig = figure()
ax = fig[:add_subplot](1,1,1)
# plot the "results" versus K (i.e., N) for each value of s_u
# i might see if this is better with JuliaPlots
#for (si, s) in enumerate(su_vals)
#    tmp_label = L"$s_u = $"*string(s)
for pars in parsets
    tmp_x = []
    tmp_y = []
    for NΓ in NΓ_vals
        if (NΓ, pars[1],pars[2],pars[3],pars[4]) in keys(NΓ_results)
            push!(tmp_x, NΓ)
            push!(tmp_y, NΓ_results[(NΓ, pars[1],pars[2],pars[3],pars[4])])
        end
        scatter(tmp_x,tmp_y)
    end
    #plot(NΓ_vals, [NΓ_results[(x,pars[1],pars[2],pars[3],pars[4])] for x in NΓ_vals])#, label=tmp_label)
end
#ax[:set_xlim](xmin=10.0^2,xmax=10.0^5.5)
#ax[:set_ylim](ymin=0.0,ymax=1.0)

#vlines([10.0^2.45,10^4.3],0,1,linestyles="dashed")
#xlabel("Population size")
#ylabel("Valley crossing probability")
#legend(loc=3)
#tight_layout()
#savefig("figs/ochs_desai_1.pdf")
