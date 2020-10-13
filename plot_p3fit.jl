# plot all relevant data 
using PyPlot, JLD2
rc("font", family="serif")
include("decompose_ttvs.jl")
@load "OUTP3/p3_fittestparams.jld2" #param_p3 lprob_p3 lprob_best pbest ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
@load "mcmc_resultstest.jld2" #par_mcmc, lprob_mcmc, nwalkers, nsteps, accept, iburn
label = "test"

pair_ttvs = decompose_ttvs(nplanet, ntrans, pbest_global)
time1 = collect(pbest_global[3] .+ range(0,stop=ntrans[1]-1,length=ntrans[1]) .* pbest_global[2])
time2 = collect(pbest_global[8] .+ range(0,stop=ntrans[2]-1,length=ntrans[2]) .* pbest_global[7])
p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
nparam=length(pbest_global)

# Make plot of decomposed simulated transit timing variations in minutes for years observed
figsize=(10,8)
subplot(211)
plot((ttmodel[1:ntrans[1]].-pbest_global[3])./365.25,pair_ttvs[1,2,1:ntrans[1]].* (24 * 60),linewidth=1.5)
plot((ttmodel[1:ntrans[1]].-pbest_global[3])./365.25,pair_ttvs[1,3,1:ntrans[1]].* (24 * 60),linewidth=1.5,color=:firebrick)
errorbar((ttmodel[1:ntrans[1]].-pbest_global[3])./365.25,(tt[1:ntrans[1]].-time1).* (24 * 60), sigtt[1:ntrans[1]].* (24 * 60),color="black",fmt=".")
ylabel("Venus TTVs (minutes)")
subplot(212)
plot((ttmodel[ntrans[1]+1:ntrans[1]+ntrans[2]].-pbest_global[8])./365.25,pair_ttvs[2,1,1:ntrans[2]].* (24 * 60),linewidth=1.5,color=:orange)
plot((ttmodel[ntrans[1]+1:ntrans[1]+ntrans[2]].-pbest_global[8])./365.25,pair_ttvs[2,3,1:ntrans[2]].* (24 * 60),linewidth=1.5,color=:firebrick)
errorbar((ttmodel[ntrans[1]+1:ntrans[1]+ntrans[2]].-pbest_global[8])./365.25,(tt[ntrans[1]+1:ntrans[1]+ntrans[2]].-time2).* (24 * 60), sigtt[ntrans[1]+1:ntrans[1]+ntrans[2]] .* (24 * 60),color="black", fmt=".")
# legend(bbox_to_anchor=(1.05, 1), loc=2,borderaxespad=0.0)
ylabel("Earth TTVs (minutes)")
xlabel("Years Observed (N)")
name = string("IMAGES/best3planetfit",label,".png")
# savefig(name)
clf()

# Make plot of best planet 3 period likelihood
figsize=(8,6)
plot(p3/365.25,exp.((lprob_p3 .-maximum(lprob_p3)))) 
# plot!( -5:8, (-5:8).^2, inset = (1, bbox(0.1,0.0,0.4,0.4)), subplot = 2)
xlabel("Period of planet 3 [years]")
ylabel("Likelihood")
name = string("IMAGES/bestp3likelihood",label,".png")
# savefig(name)
clf()

# Make plot of parameters at all MCMC steps
#   par_mcmc = zeros(nwalkers,nsteps,nparam)
#   lprob_mcmc = zeros(nwalkers,nsteps)
pname = ["mu_1","P_1","t01","e1 cos(om1)","e1 sin(om1)",
          "mu_2","P_2","t02","e2 cos(om2)","e2 sin(om2)",
          "mu_3","P_3","t03","e3 cos(om3)","e3 sin(om3)"]
function plot_MCMCsteps()
figsize=(9,5)
for i=1:5
    subplot(5,1,i)
    for j=1:nwalkers 
        plot(par_mcmc[j,1:nsteps,i])
        ylabel(pname[i])
    end
    tight_layout()
end
name = string("IMAGES/MCMCsteps",label,"p1.png")
savefig(name)
clf()
figsize=(9,5)
for i=1:5
    subplot(5,1,i)
    for j=1:nwalkers 
        plot(par_mcmc[j,1:nsteps,i+5])
        ylabel(pname[i+5])
    end
    tight_layout()
end
name = string("IMAGES/MCMCsteps",label,"p2.png")
savefig(name)
clf()
figsize=(9,5)
for i=1:5
    subplot(5,1,i)
    for j=1:nwalkers 
        plot(par_mcmc[j,1:nsteps,i+10])
        ylabel(pname[i+10])
    end
    tight_layout()
end
name = string("IMAGES/MCMCsteps",label,"p3.png")
savefig(name)
clf()
end

# Make plot of MCMC parameters after burn-in
# figsize=(8,6)
# for i=2:nparam
#   for j=1:i-1
#     scatter(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]),vec(par_mcmc[1:nwalkers,iburn:nsteps,j]))
#     xlabel(pname[i])
#     ylabel(pname[j])
#   end
# end
# name = string("IMAGES/MCMCparams",label,".png")
# savefig(name)