function show_args(args)
@show args
end
include("sim_times.jl")
include("fit_mysterybody.jl")
include("MCMC.jl")
using JLD2
###### NYR40 ########
jd1 = 2.4332825e6
jd2 = 2.4478925e6
jdsize = 1000
p3in = 4230.0
p3out = 4430.0
np3 = 100
nphase = 36
dpin = 0.0 
dpout = 2*pi
ndp = 72
nwalkers = 50
nsteps = 1000

show_args(ARGS)
label, sigma = ARGS[1],parse(Float64,ARGS[2])

function p3_run(label,sigma)
sim_times(jd1,jd2,jdsize,true,sigma,true)
datafile = string("INPUTS/tt_data",sigma,"sEMB.txt")
fit_planet3(datafile,label,jd1,jd2,jdsize,p3in,p3out,np3,nphase,true,sigma,true)
end

function moon_run(label,sigma)
sim_times(jd1,jd2,jdsize,true,sigma,false)
datafile = string("INPUTS/tt_data",sigma,"snoEMB.txt")
fit_moon(datafile,label,jd1,jd2,jdsize,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma,false)
end

p3_run(label,sigma)
moon_run(label,sigma)

function p3_mcmc(nwalkers,nsteps)
fitfile = string("OUTPUTS/p3_fit",label,"params.jld2")
foutput = string("NYR40/p3_",label)
f = jldopen(String(fitfile),"r")
MCMC(f["pbest_global"],f["lprob_best"],foutput,nsteps,nwalkers,f["nplanet"],f["ntrans"],f["tt0"],f["tt"],f["sigtt"],true,true)
end

function moon_mcmc(nwalkers,nsteps)
fitfile = string("OUTPUTS/moon_fit",label,"params.jld2")
foutput = string("NYR40/moon_",label)
f = jldopen(String(fitfile),"r")
MCMC(f["pbest_global"],f["lprob_best"],foutput,nsteps,nwalkers,f["nplanet"],f["ntrans"],f["tt0"],f["tt"],f["sigtt"],false,true)
end

p3_mcmc(nwalkers,nsteps)
moon_mcmc(nwalkers,nsteps)
