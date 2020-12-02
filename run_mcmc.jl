function show_args(args)
@show args
end
include("MCMC.jl")

# @load "OUTPUTS/p3_fittry001params.jld2" #param_p3 lprob_p3 lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase

#nwalkers = 50
#nsteps = 10000
function p3_mcmc(datafile,foutput,nwalkers,nsteps)
@load datafile lprob_best pbest_global nplanet ntrans tt0 tt sigtt 
@time lprob_mcmc,par_mcmc = MCMC(pbest_global,foutput,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,true) #no EMB yes sigsys
end

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
