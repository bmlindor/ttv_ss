# Julia v1.1
using DelimitedFiles, JLD2
include("sim_times.jl")
include("fit_mysteryplanet3.jl")
include("run_mcmc.jl")
include("fit_moon.jl")
# sig_grid = [10.0, 15.0, 30.0, 45.0, 60.0]
# years = [15, 30, 50]
# p3in = 500.0; p3out = 18000.0; np3 = 1000
# Modify the following variables as necessary:
label = "test"
jd1 = 2.4332825e6
jd2 = 2.4515445e6
jdsize = 1000
sigma = 10.0 
p3in = 4230.0
p3out = 4430.0
np3 = 10
nphase = 6
dpin = 0.0 
dpout = 2*pi
ndp = 36
nsteps = 1000
nwalkers = 50
p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
# function run_fitp3(label, p3in, p3out, np3, nphase, sigma)
sim = sim_times(jd1,jd2,jdsize,true,sigma,false)
# file = string("INPUTS/tt_data",sigma,"s.txt")
# # fit = fit_mysteryplanet3(file,label,p3in,p3out,np3,nphase,true,sigma, true)
# fit = fit_moon(file,label,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma)
# @load "OUTPUTS/p3_fittestparams.jld2" #param_p3 lprob_p3 lprob_best pbest ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
# par_mcmc, lprob_mcmc = MCMC(pbest, label, nsteps, nwalkers, nplanet, ntrans, tt0, tt, sigtt) 
@load "OUTPUTS/moon_fittestparams.jld2" #pbest_dp lprob_dp lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase dpin dpout phiphase
par_mcmc, lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,true) 

# file = string("OUTPUTS/mcmc",label,".jld2")
# @load "OUTPUTS/mcmcresultstest.jld2" #par_mcmc, lprob_mcmc, accept, iburn, steps, nwalkers, nsteps

# @load "OUTPUTS/p3_fit_test.jld2"
# param_p3 chi_p3 chi_best pbest tt0 tt ttmodel sigtt p3in p3out np3 nphase
