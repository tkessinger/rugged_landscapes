#!/usr/bin/env julia

## plt_ochs_desai.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## Open outputs of ochs_desai_sim.jl and plot results.


# tell Julia where the module is located
# include("pop_sim.jl")
using PyPlot, JLD, Glob, PopSim, LaTeXStrings

files = glob("output/ochs_desai_sims/ochs_desai_sims_*.jld")

Ns_results = Dict{Tuple{Int64, Float64}, Float64}()

for (fi, file) in enumerate(files)
    f = load(file)
    results = f["results"]
    params = f["params"]
    Ns_results[params.K,params.slist[1]] = results
end
K_vals = sort(unique(x[1] for x in keys(Ns_results)))
su_vals = sort(unique(x[2] for x in keys(Ns_results)))

fig = figure()
ax = fig[:add_subplot](1,1,1)
for (si, s) in enumerate(su_vals)
    tmp_label = L"$s_u = $"*string(s)
    semilogx(K_vals, [Ns_results[(x,s)] for x in K_vals], label=tmp_label)
end
ax[:set_xlim](xmin=10.0^2,xmax=10.0^5.5)
ax[:set_ylim](ymin=0.0,ymax=1.0)

vlines([10.0^2.45,10^4.3],0,1,linestyles="dashed")
#vlines(10.0^4.3,0,1,ls="--")
xlabel("Population size")
ylabel("Valley crossing probability")
legend(loc=3)
tight_layout()
savefig("figs/ochs_desai_1.pdf")
