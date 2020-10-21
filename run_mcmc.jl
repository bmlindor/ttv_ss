include("MCMC.jl")
# using PyPlot,Unitful,UnitfulAstro,LinearAlgebra

# @load "OUTPUTS/p3_fittestparams.jld2" #param_p3 lprob_p3 lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
# label = "testing"
@load "OUTPUTS/moon_fittestparams.jld2"
label = "systest"
nsteps = 1000
nwalkers = 50
# MCMC(param::Array{Float64,1},label::String,
#   nsteps::Int64,nwalkers::Int64,nplanet::Int64,ntrans::Array{Int64,1},
#   tt0::Array{Float64,1},tt::Array{Float64,1},sigtt::Array{Float64,1},
#   EMB::Bool=true,use_sigsys::Bool=true) we *are* using sigsys

# @time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,true,true) #yes EMB yes sigsys
# @time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,true,false) #yes EMB no sigsys
@time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,true) #no EMB yes sigsys
# @time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,false) #no EMB no sigsys
# param = [pbest;1e-4^2]
# nsteps = 1000
# nwalkers = 50
# nplanet = 3
# ntrans = [82,51,2]
# par_mcmc,chi_mcmc = MCMC(param,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt)