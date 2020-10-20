# Julia v1.1
using DelimitedFiles, JLD2
include("sim_times.jl")
include("run_mcmc.jl")
include("fit_mysteryplanet3.jl")
# sig_grid = [10.0, 15.0, 30.0, 45.0, 60.0]
# years = [15, 30, 50]
# p3in = 500.0; p3out = 18000.0; np3 = 1000
# Modify the following variables as necessary:
label = "try000"
jd1 = 2.4332825e6
jd2 = 2.4515445e6
jdsize = 1000
sigma = 30.0 
p3in = 4230.0
p3out = 4430.0
np3 = 10
nphase = 6
nsteps = 100000
nwalkers = 50
p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)

@time sim = sim_times(jd1,jd2,jdsize,true,sigma)
file = string("INPUTS/tt_data",sigma,"snoEMB.txt")
@time fit = fit_mysteryplanet3(file,label,p3in,p3out,np3,nphase,true,sigma, true)
@load "OUTP3/p3_fittestparams.jld2" #param_p3 lprob_p3 lprob_best pbest ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
@time par_mcmc, lprob_mcmc = MCMC(pbest_global, label, nsteps, nwalkers, nplanet, ntrans, tt0, tt, sigtt) #par_mcmc, lprob_mcmc, accept, iburn, steps, nwalkers, nsteps

