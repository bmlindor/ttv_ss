# Julia v1.1
include("sim_times.jl")
include("fit_mysterybody.jl")
include("MCMC.jl")

jd1 = 2.4332825e6
jd2 = 2.4515445e6
jdsize = 1000
# p3in = 4230.0
# p3out = 4430.0
# p3in = 1000.0
sigma = 30.0
# p3in = 500.0
# p3out = 5000.0
np3 = 100
nphase = 100
dpin = 0.0 
dpout = 2*pi
ndp = 72

nsteps = 10000
nwalkers = 50

function full_run()
# Modify the following variables as necessary:
label = "try01"

function full_run(label,sigma)
sim_times(jd1,jd2,jdsize,true,sigma,true)
datafile = string("INPUTS/tt_data",sigma,"sEMB.txt")
fit_planet3(datafile,label,jd1,jd2,jdsize,p3in,p3out,np3,nphase,true,sigma,true)
@load ("OUTPUTS/p3_fittry01params.jld2")
# @load "OUTP3/p3_fittestparams.jld2" #param_p3 lprob_p3 lprob_best pbest ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
@time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt) #par_mcmc,lprob_mcmc,accept,iburn,steps,nwalkers,nsteps
end

function full_moonrun()
# label = ["mtry1","mtry2","mtry3","mtry4","mtry5","mtry6","mtry7"]
# sigma = [10.0, 15.0, 30.0, 45.0, 60.0, 120.0, 240.0]

for i=1:length(sigma)
	@time sim = sim_times(jd1,jd2,jdsize,true,sigma[i],false)
	file = string("INPUTS/tt_data",sigma[i],"snoEMB.txt")
	@time lprobfit,bestfit = fit_moon(file,label[i],jd1,jd2,jdsize,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma[i],false)
end
# @load "OUTPUTS/moon_fittestparams.jld2" #pbest_dp lprob_dp lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase dpin dpout phiphase
# @time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,true) 
end
