function show_args(args)
@show args
end
using Profile
using DelimitedFiles,JLD2,LsqFit,Statistics
if !@isdefined(TTVFaster)
    include("TTVFaster/src/TTVFaster.jl")
  using Main.TTVFaster
end
import Main.TTVFaster.ttv_wrapper
import Main.TTVFaster.chisquare
include("regress.jl")
include("sim_times.jl")
include("fit_planet2.jl")
include("fit_planet3.jl")
include("fit_planet4.jl")
include("fit_planet5.jl")
include("fit_moon.jl")
include("MCMC.jl")

show_args(ARGS)
label,sigma,nyear,runtype,obs,nsteps=ARGS[1],parse(Float64,ARGS[2]),parse(Float64,ARGS[3]),ARGS[4],ARGS[5],parse(Int64,ARGS[6])
function parse_model(label::String)
	nplanet = 0
	nmoon = 0
	for obj in 1:length(label)
	  if label[obj]=='p'
	    # println(obj)
	    nplanet += 1
	  end
	  if label[obj]=='m'
	  	nmoon += 1
	  end
	end
	return nplanet,nmoon
end
nplanet,nmoon=parse_model(label)
# println(nplanet," planets and ", nmoon," moons in model")
# Initialize variables and period ranges
jd1=2.4332825e6
tref=2430000
tol=1e-5
nphase=36 #wide: 100,36,180 
p3in,p3out,np3=10*365.25,15*365.25,100
dpin,dpout,ndp=2.1,2.4,100
p4in,p4out,np4=1.8*365.25,2.2*365.25,100
p5in,p5out,np5=28*365.25,33*365.25,100
nwalkers=75
# Change search grids to accomodate for wider distributions when time spans are shorter
# if nyear>36
# 	p3in,p3out=11.4*365.25,12.1*365.25
# 	dpin,dpout=2.25,2.37
# 	nsteps=10000
# elseif nyear>=30 && nyear<=36
# 	p3in,p3out=10.8*365.25,12.3*365.25
# 	dpin,dpout=2.2,2.38
# 	nsteps=50000
# elseif nyear>=12 && nyear<=29
# 	p3in,p3out=10.6*365.25,14.2*365.25 
# 	dpin,dpout=2.1,2.52
# 	nsteps=50000
# elseif nyear<12
# 	p3in,p3out=10*365.25,15*365.25 
# 	nsteps=100000
# end
# Run markov chains
function planet_mcmc(nplanet,nsteps,obs::String)
	if obs=="fromEMB"
	fitfile=string("FITS/fromEMB/p",nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
	foutput=string("MCMC/fromEMB/p",nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
	elseif obs=="fromEV"
	fitfile=string("FITS/p",nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
	foutput=string("MCMC/p",nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
	end
	p=jldopen(String(fitfile),"r")
	MCMC(foutput,p[string("best_p",nplanet)],p[string("lprob_best_p",nplanet)],nsteps,nwalkers,nplanet,p["ntrans"][1:nplanet],p["tt0"],p["tt"],p["sigtt"],true,true)
end 
function moon_mcmc(nplanet,nsteps,label::String)
	if label=="Hppmpp" # load p4 search after moon fit
	fitfile=string("FITS/moonp4_fit",sigma,"s",nyear,"yrs.jld2")
	foutput=string("MCMC/moonp4_mcmc",sigma,"s",nyear,"yrs.jld2")
	m=jldopen(String(fitfile),"r")
	MCMC(foutput,m["best_p4"],m["lprob_best_p4"],nsteps,nwalkers,4,m["ntrans"][1:4],m["tt0"],m["tt"],m["sigtt"],true,false)
	else
	fitfile=string("FITS/p",nplanet,"moon_fit",sigma,"s",nyear,"yrs.jld2")
	foutput=string("MCMC/p",nplanet,"moon_mcmc",sigma,"s",nyear,"yrs.jld2")
	m=jldopen(String(fitfile),"r")
	MCMC(foutput,m["best_dp"],m["lprob_best_dp"],nsteps,nwalkers,nplanet,m["ntrans"][1:nplanet],m["tt0"],m["tt"],m["sigtt"],true,false)	
	end
end
# Perform grid searches
function run_grid(label::String,obs::String)
	if label=="Hpp"
		sim_times(jd1,sigma,nyear,obs)
		fit_planet2(jd1,sigma,nyear,tref,tol,obs)
	elseif label=="Hppp"
		fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,obs)
	elseif label=="Hpppp"
		fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,obs)
	elseif label=="Hppppp"
		fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,obs)
	end
	if label=="Hppm"
		fit_planet2(jd1,sigma,nyear,tref,tol,obs)
		fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,2)
	elseif label=="Hppmp"
		fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,obs)
		fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,3)
	elseif label=="Hppmpp"
		# fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,obs)
		# fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,4)
		fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase)
	elseif label=="Hppmppp"
		fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,obs)
		fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,5)
	end
end
function test_fit(label::String,obs::String)
	nphase=36 #wide: 100,36,180 
	p3in,p3out,np3=11.4*365.25,12.2*365.25,10
	dpin,dpout,ndp=2.29,2.35,10
	p4in,p4out,np4=1.6*365.25,3*365.25,20
	p5in,p5out,np5=28*365.25,30*365.25,10
	#sim_times(jd1,sigma,nyear,obs)
	#@time fit_planet2(jd1,sigma,nyear,tref,tol,obs)
	#@time fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,obs)
	@time fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,obs)
	#@time fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,obs)
	if label=="Hppmpp" 
		@time fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,3)
		@time fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase)
	end
end
if runtype=="fittest" 
	println("Testing")
	test_fit(label,obs)
	# @time planet_mcmc(nplanet,1000,obs)
elseif	runtype=="test" && label=="Hppmpp"
	println("Testing")
	test_fit(label,obs)
	@time moon_mcmc(nplanet,1000,label)
elseif runtype=="mctest"
	println("Testing")
	# test_fit(label,obs)
	@time planet_mcmc(nplanet,1000,obs)
end
if runtype=="full" && nmoon==0
	@time run_grid(label,obs)
	@time planet_mcmc(nplanet,nsteps,obs)
elseif runtype=="full" && nmoon>0
	@time run_grid(label,obs)
	@time moon_mcmc(nplanet,nsteps,label)
end
if runtype=="grid"
	@time run_grid(label,obs)
end
if runtype=="mcmc" && nmoon==0
	@time planet_mcmc(nplanet,nsteps,obs)
elseif runtype=="mcmc" && nmoon>0
	@time moon_mcmc(nplanet,nsteps,label)
end
if runtype=="wide"
	nphase=36 #wide: 100,36,180 
	p3in,p3out,np3=5*365.25,22*365.25,100
	dpin,dpout,ndp=0.0,2*pi,100
	p4in,p4out,np4=1.5*365.25,5*365.25,100
	p5in,p5out,np5=22*365.25,36*365.25,100
	fit_wide(jd1,sigma,nyear,tref,tol,5*365.25,20*365.25,100,36,0.0,2*pi,180,1.5*365.25,5*365.25,100,obs)
	# @time fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,obs)
	# @time fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,obs)
	# @time fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,obs)
	# @time fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,3)
	# planet_mcmc(3,obs)
	# moon_mcmc(3,label)
end

# elseif runtype=="wide" && label=="test"
# 	wide_run(true,36,3*365.25,30*365.25,200,0.0,2pi,180)
# 	println("No routine available with that runtype and/or label.")

# function full_moonrun()
# # label=["mtry1","mtry2","mtry3","mtry4","mtry5","mtry6","mtry7"]
# # sigma=[10.0, 15.0, 30.0, 45.0, 60.0, 120.0, 240.0]
# for i=1:length(sigma)
# 	@time sim=sim_times(jd1,jd2,jdsize,true,sigma[i],false)
# 	file=string("INPUTS/tt_data",sigma[i],"snoEMB.txt")
# 	@time lprobfit,bestfit=fit_moon(file,label[i],jd1,jd2,jdsize,p3in,p3out,np3,nphase,dpin,dpout,ndp,true,sigma[i],false)
# end
# end
