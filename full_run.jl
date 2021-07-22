function show_args(args)
@show args
end
using Profile
include("sim_times.jl")
include("fit_planet3.jl")
include("fit_moon.jl")
include("MCMC.jl")
show_args(ARGS)
runtype, label = ARGS[1], ARGS[4]
sigma, nyear = parse(Float64,ARGS[2]),parse(Float64,ARGS[3])
# Initialize variables and fixed values
jd1 = 2.4332825e6
np3,nphase,ndp = 200,36,72 #100,36,180 
p3in,p3out,dpin,dpout = zeros(4)
nwalkers,nsteps = 50,0
# Change search grids to accomodate for wider distributions when time spans are shorter
if nyear==40
	p3in,p3out=11.4*365.25,12.1*365.25
	dpin,dpout=2.25,2.37
	nsteps=10000
end
if nyear==30
	p3in,p3out=10.8*365.25,12.3*365.25
	dpin,dpout=2.2,2.38
	nsteps=50000
end
if nyear==20 || nyear==25
	p3in,p3out=10.6*365.25,14.2*365.25 
	dpin,dpout=2.1,2.52
	nsteps=50000
end
if nyear==15
	p3in,p3out=10.6*365.25,14.2*365.25 
	dpin,dpout=2.1,2.52
	nsteps=50000
end
if nyear==10
	p3in,p3out=10*365.25,15*365.25 
	nsteps=50000
end
# Planet 3 detection and characterization routine
if label=="ppp"
	sim_times(jd1,nyear,true,sigma,true) 
	function grid_run(p3in,p3out,np3,nphase)
		datafile = string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt")
		# fileout = string("FITS/p3_wide")
		@time fit_planet3(datafile,jd1,nyear,p3in,p3out,np3,nphase,true,sigma,true)
	end
	function run_mcmc()
		fitfile = string("FITS/p3_fit",sigma,"s",nyear,"yrs.jld2")
		foutput = string("MCMC/p3_mcmc",sigma,"s",nyear,"yrs.jld2")
		p = jldopen(String(fitfile),"r")
		@time MCMC(foutput,p["pbest_global"],p["lprob_best"],nsteps,nwalkers,p["nplanet"],p["ntrans"],p["tt0"],p["tt"],p["sigtt"],true,true)	
	end    
	if runtype=="grid"
		grid_run(p3in,p3out,np3,nphase)
	elseif runtype=="mcmc"
		run_mcmc()
	elseif runtype=="full"
		grid_run(p3in,p3out,np3,nphase)
		run_mcmc()
	elseif runtype=="wide"
		grid_run(1000.0,5500.0,100,36)
	end
end

# Moon detection and characterization routine
if label=="ppmp"
	  sim_times(jd1,nyear,true,sigma,false)
	  function grid_run(p3in,p3out,np3,nphase,dpin,dpout,ndp)
		  datafile = string("INPUTS/tt_",sigma,"snoEMB",nyear,"yrs.txt")
			@time fit_moon(datafile,jd1,nyear,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma,false)   
	  end
	  function run_mcmc()
		  fitfile = string("FITS/moon_fit",sigma,"s",nyear,"yrs.jld2")
			foutput = string("MCMC/moon_mcmc",sigma,"s",nyear,"yrs.jld2")
			m = jldopen(String(fitfile),"r")
			@time MCMC(foutput,m["pbest_global"],m["lprob_best"],nsteps,nwalkers,m["nplanet"],m["ntrans"],m["tt0"],m["tt"],m["sigtt"],true,false)	  
	  end
  if runtype=="grid"
  	grid_run(p3in,p3out,np3,nphase,dpin,dpout,ndp)
	elseif runtype=="mcmc"  
   	run_mcmc()
  elseif runtype=="full"
   	grid_run(p3in,p3out,np3,nphase,dpin,dpout,ndp)
   	run_mcmc()
  elseif runtype=="wide"
		grid_run(1000.0,5000.0,100,36,0.0,2pi,180)
  end  
end

# if runtype=="sim"
# 	if label=="ppp"
# 	  sim_times(jd1,nyear,true,sigma,true) 
# 	end
# 	if label=="ppmp"
# 	  sim_times(jd1,nyear,true,sigma,false)
# 	end
# end
# if runtype=="grid" 
# 	if label=="ppp"
# 		# sim_times(jd1,nyear,true,sigma,true)
# 		datafile = string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt")
# 		@time fit_planet3(datafile,jd1,nyear,p3in,p3out,np3,nphase,true,sigma,true)
# 	end
# 	if label=="ppmp"
#     # sim_times(jd1,nyear,true,sigma,false)
# 		datafile = string("INPUTS/tt_",sigma,"snoEMB",nyear,"yrs.txt")
# 		@time fit_moon(datafile,jd1,nyear,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma,false)
# 	end
# end
# if runtype=="mcmc" 
# 	if label=="ppp"
# 		fitfile = string("FITS/p3_fit",sigma,"s",nyear,"yrs.jld2")
# 		foutput = string("MCMC/p3_mcmc",sigma,"s",nyear,"yrs.jld2")
# 		p = jldopen(String(fitfile),"r")
# 		MCMC(foutput,p["pbest_global"],p["lprob_best"],nsteps,nwalkers,p["nplanet"],p["ntrans"],p["tt0"],p["tt"],p["sigtt"],true,true)
# 	end
# 	if label=="ppmp"
#     fitfile = string("FITS/moon_fit",sigma,"s",nyear,"yrs.jld2")
# 		foutput = string("MCMC/moon_mcmc",sigma,"s",nyear,"yrs.jld2")
# 		m = jldopen(String(fitfile),"r")
# 		@time MCMC(foutput,m["pbest_global"],m["lprob_best"],nsteps,nwalkers,m["nplanet"],m["ntrans"],m["tt0"],m["tt"],m["sigtt"],true,false)
# 	end 
# end
# if runtype=="grid" && label=="ppmp"
# sim_times(jd1,nyear,true,sigma,false)
# datafile = string("INPUTS/tt_",sigma,"snoEMB",nyear,"yrs.txt")
# @time fit_moon(datafile,jd1,nyear,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma,false)
# end
# if runtype=="mcmc" && label=="ppmp"
# fitfile = string("FITS/moon_fit",sigma,"s",nyear,"yrs.jld2")
# foutput = string("MCMC/moon_mcmc",sigma,"s",nyear,"yrs.jld2")
# f = jldopen(String(fitfile),"r")
# @time MCMC(foutput,f["pbest_global"],f["lprob_best"],nsteps,nwalkers,f["nplanet"],f["ntrans"],f["tt0"],f["tt"],f["sigtt"],true,false)
# end

# function p4_run(sigma,nyear)
# if runtype=="ppp"
# sim_times(jd1,nyear,jdsize,true,sigma,true)
# datafile = string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt")
# @time fit_planet4(datafile,jd1,nyear,jdsize,p3in,p3out,np3,nphase,p4in,p4out,np4,true,sigma,true)
# fitfile = string("FITS/p4_fit",sigma,"s",nyear,"yrs.jld2")
# foutput = string("MCMC/p4_mcmc",sigma,"s",nyear,"yrs.jld2")
# f = jldopen(String(fitfile),"r")
# MCMC(foutput,f["pbest_global"],f["lprob_best"],nsteps,nwalkers,f["nplanet"],f["ntrans"],f["tt0"],f["tt"],f["sigtt"],true,true)
# end
# end


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
