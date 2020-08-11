# Julia v1.1
using JLD2, DelimitedFiles
include("sim_times.jl")
include("fit_mysteryplanet3.jl")
# sig_grid = [10.0, 15.0, 30.0, 45.0, 60.0]
# p3in = 500.0; p3out = 18000.0; np3 = 1000
label = "test"
p3in = 4230.0
p3out = 4430.0
np3 = 10
nphase = 6
sigma = 30.0
p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
# function run_fitp3(label, p3in, p3out, np3, nphase, sigma)
	# sim = sim_times(2.4332825e6, 2.4515445e6, 1000, true, sigma)
file = string("INPUTS/tt_data",sigma,"s.txt")
fit = fit_mysteryplanet3(file,label,p3in,p3out,np3,nphase,true, sigma)
	# @load string(“OUTPUTS/p3_fit_params”,label,”.jld2”)
	# par_mcmc, lprob_mcmc = MCMC(param,nsteps,nwalkers,nplanet,ntrans,tt0, tt, sigtt) 
# end

# run_fitp3("test", p3in, p3out, np3, nphase, 30.0)

# @load "OUTPUTS/p3_fit_params.jld2"
# param = [pbest;1e-4^2]
nsteps = 1000
nwalkers = 50
nplanet = 3
ntrans = [82,51,2]
# par_mcmc, lprob_mcmc = MCMC(param,nsteps,nwalkers,nplanet,ntrans,tt0, tt, sigtt) 

# @load "OUTPUTS/p3_fit_test.jld2"
# param_p3 chi_p3 chi_best pbest tt0 tt ttmodel sigtt p3in p3out np3 nphase