function show_args(args)
@show args
end
using Profile
include("sim_times.jl")
# include("fit_planet3.jl")
# include("fit_planet4.jl")
include("fit_moon.jl")
include("MCMC.jl")
show_args(ARGS)
runtype, label = ARGS[1], ARGS[4]
sigma, nyear = parse(Float64,ARGS[2]),parse(Float64,ARGS[3])
# Initialize variables and fixed values
jd1 = 2.4332825e6
np3,nphase,ndp = 200,36,180 #wide: 100,36,180 
p3in,p3out,dpin,dpout = zeros(4)
p4in,p4out,np4= 1.5*365.25,5*365.25,10
nwalkers,nsteps = 50,1000
# Change search grids to accomodate for wider distributions when time spans are shorter
if nyear==40
	p3in,p3out=11.4*365.25,12.1*365.25
	dpin,dpout=2.25,2.37
	nsteps=10000
end
if nyear>=30 && nyear<=39
	p3in,p3out=10.8*365.25,12.3*365.25
	dpin,dpout=2.2,2.38
	nsteps=50000
end
if nyear>=20 && nyear<=29
	p3in,p3out=10.6*365.25,14.2*365.25 
	dpin,dpout=2.1,2.52
	nsteps=50000
end
if nyear>=15 && nyear<=19
	p3in,p3out=10.6*365.25,14.2*365.25 
	dpin,dpout=2.1,2.52
	nsteps=50000
end
if nyear<=15
	p3in,p3out=10*365.25,15*365.25 
	nsteps=100000
end
	# For planet 3 search routine with EMB simulations, use
	# datafile = string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt") and EMB==true
	# @time fit_planet3(datafile,jd1,nyear,p3in,p3out,np3,nphase,true,sigma,true)

# sim_times(jd1,nyear,true,sigma,false)
# Fit V + E transit times and perform 3-planet grid search
function grid_run(p3in,p3out,np3,nphase,dpin,dpout,ndp)
  datafile = string("INPUTS/tt_",sigma,"snoEMB",nyear,"yrs.txt")
	@time fit_moon(datafile,jd1,nyear,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma,false)   
end
fitfile = string("FITS/moon_fit",sigma,"s",nyear,"yrs.jld2")
# Run 3-planet markov chain
function p3_mcmc()
foutput = string("MCMC/p3_mcmc",sigma,"s",nyear,"yrs.jld2")
p = jldopen(String(fitfile),"r")
@time MCMC(foutput,p["best_p3"],p["lprob_best_p3"],nsteps,nwalkers,3,p["ntrans"][1:3],p["tt0"],p["tt"],p["sigtt"],true,true)	
end 
# Run lunar markov chain
function moon_mcmc()
	foutput = string("MCMC/moon_mcmc",sigma,"s",nyear,"yrs.jld2")
	m = jldopen(String(fitfile),"r")
	@time MCMC(foutput,m["best_dp"],m["lprob_best_dp"],nsteps,nwalkers,3,m["ntrans"][1:3],m["tt0"],m["tt"],m["sigtt"],true,false)	  
end

if runtype=="mcmc" && label=="ppmp"
	moon_mcmc()
elseif runtype=="mcmc" && label=="ppp"
	p3_mcmc()
elseif runtype=="grid" && label=="ppmp"
 	grid_run(p3in,p3out,np3,nphase,dpin,dpout,ndp)
elseif runtype=="grid" && label=="pppp" # Planet 4 detection routine
	@time fit_planet4(sigma,nyear,p4in,p4out,np4)
elseif runtype=="full" && label=="ppmp" 
 	grid_run(p3in,p3out,np3,nphase,dpin,dpout,ndp)
 	moon_mcmc()
else
	println("No routine available with that runtype or label.")
  # elseif runtype=="wide"
		# grid_run(1000.0,5000.0,100,36,0.0,2pi,180)
end

# function full_moonrun()
# # label = ["mtry1","mtry2","mtry3","mtry4","mtry5","mtry6","mtry7"]
# # sigma = [10.0, 15.0, 30.0, 45.0, 60.0, 120.0, 240.0]

# for i=1:length(sigma)
# 	@time sim = sim_times(jd1,jd2,jdsize,true,sigma[i],false)
# 	file = string("INPUTS/tt_data",sigma[i],"snoEMB.txt")
# 	@time lprobfit,bestfit = fit_moon(file,label[i],jd1,jd2,jdsize,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma[i],false)
# end
# # @load "OUTPUTS/moon_fittestparams.jld2" #pbest_dp lprob_dp lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase dpin dpout phiphase
# # @time par_mcmc,lprob_mcmc = MCMC(pbest_global,label,nsteps,nwalkers,nplanet,ntrans,tt0,tt,sigtt,false,true) 
# end
