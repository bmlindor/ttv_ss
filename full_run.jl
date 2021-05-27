function show_args(args)
@show args
end
using Profile
include("sim_times.jl")
include("fit_planet3.jl")
include("fit_moon.jl")
include("MCMC.jl")
nwalkers = 50
nsteps = 10000 #10000
jd1 = 2.4332825e6
p3in = 4163.8
p3out = 4419.5
np3 = 200 #100
nphase = 36 #36
dpin = 2.28 #0.0
dpout = 2.34 #2pi
ndp = 72 #180 #72
# p4in = 
# p4out = 
# np4 = 100
show_args(ARGS)
runtype, label = ARGS[1], ARGS[4]
sigma, nyear = parse(Float64,ARGS[2]),parse(Float64,ARGS[3])

if label=="ppp"
	sim_times(jd1,nyear,true,sigma,true) 
	function grid_run()
		datafile = string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt")
		@time fit_planet3(datafile,jd1,nyear,p3in,p3out,np3,nphase,true,sigma,true)
	end
	function run_mcmc()
		fitfile = string("FITS/p3_fit",sigma,"s",nyear,"yrs.jld2")
		foutput = string("MCMC/p3_mcmc",sigma,"s",nyear,"yrs.jld2")
		p = jldopen(String(fitfile),"r")
		@time MCMC(foutput,p["pbest_global"],p["lprob_best"],nsteps,nwalkers,p["nplanet"],p["ntrans"],p["tt0"],p["tt"],p["sigtt"],true,true)	
	end    
	if runtype=="grid"
		grid_run()
	elseif runtype=="mcmc"
		run_mcmc()
	elseif runtype=="full"
		grid_run()
		run_mcmc()
	end
    
end

if label=="ppmp"
	  sim_times(jd1,nyear,true,sigma,false)
	  function grid_run()
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
  	grid_run()
	elseif runtype=="mcmc"  
   	run_mcmc()
   elseif runtype=="full"
   	grid_run()
   	run_mcmc()
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
