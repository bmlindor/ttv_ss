include("MCMC.jl")
# using PyPlot,Unitful,UnitfulAstro,LinearAlgebra

# @load "OUTPUTS/p3_fittestparams.jld2" #param_p3 lprob_p3 lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
# label = "testing"
@load ("OUTPUTS/moon_fittestparams.jld2")
label = "mtest"
nsteps = 1000
nwalkers = 50

# @time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,true,true) #yes EMB yes sigsys
@time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,true) #no EMB yes sigsys
# param = [pbest;1e-4^2]
# nsteps = 1000
# nwalkers = 50
# nplanet = 3
# ntrans = [82,51,2]
# par_mcmc,chi_mcmc = MCMC(param,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt)