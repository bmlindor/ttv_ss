# Julia v1.1
using DelimitedFiles,JLD2
include("sim_times.jl")
# include("mcmc.jl")
include("fit_moon.jl")
# sig_grid = [10.0,15.0,30.0,45.0,60.0]
# years = [15,30,50]
# p3in = 500.0; p3out = 18000.0; np3 = 1000
# Modify the following variables as necessary:
label = "mtry1"
jd1 = 2.4332825e6
jd2 = 2.4515445e6
jdsize = 1000
sigma = 30.0 
p3in = 4230.0
p3out = 4430.0
np3 = 20
nphase = 10
dpin = 0.0 
dpout = 2*pi
ndp = 36
nsteps = 3000
nwalkers = 50

# @time sim = sim_times(jd1,jd2,jdsize,true,sigma,false)
file = string("INPUTS/tt_data",sigma,"snoEMB.txt")
@time fit = fit_moon(file,label,jd1,jd2,jdsize,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma,false)
# @load "OUTPUTS/moon_fittestparams.jld2" #pbest_dp lprob_dp lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase dpin dpout phiphase
# @time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,true) 