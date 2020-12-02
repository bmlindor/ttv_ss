function show_args(args)
@show args
end
include("MCMC.jl")

<<<<<<< HEAD
# @load "OUTPUTS/p3_fittry001params.jld2" #param_p3 lprob_p3 lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase

#nwalkers = 50
#nsteps = 10000
function p3_mcmc(datafile,foutput,nwalkers,nsteps)
@load datafile lprob_best pbest_global nplanet ntrans tt0 tt sigtt 
@time lprob_mcmc,par_mcmc = MCMC(pbest_global,foutput,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,true) #no EMB yes sigsys
end
=======
# file1 = "OUTPUTS/p3_fittry001params.jld2"
# extractdata(file1)
@load "OUTPUTS/p3_fittry001params.jld2" #param_p3 lprob_p3 lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
label = "try001"
nsteps = 100000
nwalkers = 50
@time lprob_mcmc,par_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,true,true) #yes EMB yes sigsys

#@load ("OUTPUTS/moon_fitmtestparams.jld2")
#label = "mtry1"
#nsteps = 200000
#nwalkers = 50
#@time lprob_mcmc,par_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,true) #no EMB yes sigsys
>>>>>>> 3ce36cb5e6f823a4ad9d2e3ac1c1fbb614f895a3

function moon_mcmc(datafile,foutput,nwalkers,nsteps)
@load datafile lprob_best pbest_global nplanet ntrans tt0 tt sigtt
@time lprob_mcmc,par_mcmc = MCMC(pbest_global,foutput,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,true,true) #yes EMB yes sigsys
end
# param = [pbest;1e-4^2]
# nsteps = 1000
# nwalkers = 50
# nplanet = 3
# ntrans = [82,51,2]
# par_mcmc,chi_mcmc = MCMC(param,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt)
