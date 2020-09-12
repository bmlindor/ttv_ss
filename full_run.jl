# Julia v1.1
using DelimitedFiles, JLD2
include("sim_times.jl")
include("fit_mysteryplanet3.jl")
include("run_mcmc.jl")
include("fit_moon.jl")
# sig_grid = [10.0, 15.0, 30.0, 45.0, 60.0]
# p3in = 500.0; p3out = 18000.0; np3 = 1000
# Modify the following variables as necessary:
label = "test"
p3in = 4230.0
p3out = 4430.0
np3 = 10
nphase = 6
dpin = 0.0 # --> 0
dpout = 2*pi
ndp = 6 #36
nphiphase = 6
sigma = 30.0 
nsteps = 1000
nwalkers = 50
p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
# function run_fitp3(label, p3in, p3out, np3, nphase, sigma)
sim = sim_times(2.4332825e6, 2.4515445e6, 1000, true, sigma, false)
file = string("INPUTS/tt_data",sigma,"s.txt")
# fit = fit_mysteryplanet3(file,label,p3in,p3out,np3,nphase,true,sigma, true)
fit = fit_moon(file,label,p3in,p3out,np3,nphase,dpin, dpout, ndp, nphiphase,true,sigma)
# @load "OUTPUTS/p3_fittestparams.jld2" #param_p3 lprob_p3 lprob_best pbest ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
# par_mcmc, lprob_mcmc = MCMC(pbest, label, nsteps, nwalkers, nplanet, ntrans, tt0, tt, sigtt) 
# file = string("OUTPUTS/mcmc",label,".jld2")
# @load "OUTPUTS/mcmcresultstest.jld2" #par_mcmc, lprob_mcmc, accept, iburn, steps, nwalkers, nsteps
# par_mcmc, lprob_mcmc = MCMC(param,nsteps,nwalkers,nplanet,ntrans,tt0, tt, sigtt) 

# @load "OUTPUTS/p3_fit_test.jld2"
# param_p3 chi_p3 chi_best pbest tt0 tt ttmodel sigtt p3in p3out np3 nphase