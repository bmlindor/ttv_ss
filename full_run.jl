function show_args(args)
@show args
end
using Profile,TTVFaster
using DelimitedFiles,JLD2,DataFrames,CSV,LsqFit,Statistics
include("regress.jl")
include("sim_times.jl")
include("fit_planet2.jl")
include("fit_planet3.jl")
include("fit_planet4.jl")
include("fit_planet5.jl")
include("fit_moon.jl")
include("MCMC.jl")

show_args(ARGS)
label,runtype,obs,nper=ARGS[1],ARGS[2],ARGS[3],parse(Int64,ARGS[4])
#parse(Float64,ARGS[2]),parse(Float64,ARGS[3]),ARGS[4],ARGS[5],parse(Int64,ARGS[6])
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
tref=2430000; tol=1e-5
nphase=36 #wide: 100,36,180 
np3,np4,np5=[nper, nper, nper]

# Run markov chains
#sigma=10; nyear=20
nwalkers=75;nsteps=50000
function planet_mcmc(sigma,nyear,nplanet,nsteps,obs::String)
	if obs=="fromEMB"
	fitfile=string("FITS/fromEMB/p",nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
	foutput=string("MCMC/fromEMB/p",nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
	results=string("MCMC/fromEMB/p",nplanet,"_results",sigma,"s",nyear,"yrs.txt")
	elseif obs=="fromEV"
	fitfile=string("FITS/p",nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
	foutput=string("MCMC/p",nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
	results=string("MCMC/p",nplanet,"_results",sigma,"s",nyear,"yrs.txt")
	end
	p=jldopen(String(fitfile),"r")
	MCMC(foutput,p[string("best_p",nplanet)],p[string("lprob_best_p",nplanet)],nsteps,nwalkers,nplanet,p["ntrans"][1:nplanet],p["tt0"],p["tt"],p["sigtt"],true,true,results)
end 
function moon_mcmc(sigma,nyear,nplanet,nsteps,label::String)
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
# Change search grids to accomodate for wider distributions when time spans are shorter

function run_grid(sigma,nyear,label::String,obs::String)
 if nyear>=26
 	p3in,p3out=11.2*365.25,11.9*365.25
	p4in,p4out=1.8*365.25,1.95*365.25
 	dpin,dpout=2.25,2.37
 elseif 22<=nyear<26
 	p3in,p3out=11.2*365.25,12*365.25
	p4in,p4out=1.8*365.25,2*365.25
 	dpin,dpout=2.2,2.38
 elseif 18<=nyear<22
 	p3in,p3out=11.2*365.25,12.4*365.25
	p4in,p4out=1.8*365.25,2.2*365.25
 	dpin,dpout=2.2,2.38
 elseif 12<=nyear<18
 	p3in,p3out=11*365.25,12.8*365.25 
	p4in,p4out=1.8*365.25,3.6*365.25
 	dpin,dpout=2.1,2.52
 	nsteps=50000
 elseif nyear<12
 	p3in,p3out=11*365.25,13*365.25 
	p4in,p4out=1.6*365.25,3.8*365.25
 end

	if label=="Hpp" 
	#sim_times(jd1,sigma,nyear,obs)
 	fit_planet2(jd1,sigma,nyear,tref,tol,[obs],true)
	elseif label=="Hppp"
	fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,[obs,"p3"],true)
	elseif label=="Hpppp"
  fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,[obs,"p4"],true)
	elseif label=="Hppppp"
  fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,[obs,"p5"],true)
	end
	if label=="Hppm"
 	fit_planet2(jd1,sigma,nyear,tref,tol,[obs],true)
	fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,2,"p2moon")
	elseif label=="Hppmp"
	fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,[obs,"p3"],true)
	fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,3,"p3moon",true)
	elseif label=="Hppmpp"
  #fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,[obs,"p4"],true)
	# fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,4,"p4moon",true)
	fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,"p3moonp4",true)
	elseif label=="Hppmppp"
  fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,[obs,"p5"],true)
	fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,5,"p5moon",true)
	end
end
function test_fit(sigma,nyear,label::String,obs::String)
	nphase=36 #wide: 100,36,180 
	p3in,p3out,np3=11.4*365.25,12.2*365.25,20
	dpin,dpout,ndp=2.29,2.35,10
	p4in,p4out,np4=1.6*365.25,3*365.25,20
	p5in,p5out,np5=28*365.25,30*365.25,20
	#jd1=2444239.5 ; tref =2444000
	if label=="Hpp" 
	sim_times(jd1,sigma,nyear,obs)
 	best_p2,lprob_best_p2=fit_planet2(jd1,sigma,nyear,tref,tol,[obs],true)
	elseif label=="Hppp"
	fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,[obs,"p3"],true)
	elseif label=="Hpppp"
  fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,[obs,"p4"])
	elseif label=="Hppppp"
  fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,[obs,"p5"])
	end
	if label=="Hppmpp" 
	fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,3,"p3moon")
	fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,"p3moonp4")
	end
end

if runtype=="wide"
	nphase=36 #wide: 100,36,180 
	p3in,p3out,np3=5*365.25,22*365.25,500
	dpin,dpout,ndp=0.0,2*pi,100
	p4in,p4out,np4=1.5*365.25,5*365.25,200
	p5in,p5out,np5=22*365.25,36*365.25,200
  fit_planet2(jd1,sigma,nyear,tref,tol,[obs])
  fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,[obs,"widep3"],true)
  fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,[obs,"widep4"],true)
  fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,[obs,"widep5"],true)
	# @time fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,3)
end
nyears=[15,17,19,21,23,25,27,29]#,14,13,12,11,10]16,18,20,22,24,26,28,30]#
sigmas=[10,30,60,80,90,100,110,120]
for sig in sigmas
for yr in nyears
if runtype=="full" && nmoon==0
	@time run_grid(sig,yr,label,obs)
	@time planet_mcmc(sig,yr,nplanet,nsteps,obs)
elseif runtype=="full" && nmoon>0
	@time run_grid(sig,yr,label,obs)
	@time moon_mcmc(sig,yr,nplanet,nsteps,label)
end
if runtype=="mcmc" && nmoon==0
	@time planet_mcmc(sig,yr,nplanet,nsteps,obs)
elseif runtype=="mcmc" && nmoon>0
	@time moon_mcmc(sig,yr,nplanet,nsteps,label)
end
if runtype=="grid"
run_grid(sig,yr,label,obs)
end
end
end
if runtype=="test"
nyears=[30]
sigmas=[10]
for sig in sigmas
for yr in nyears
test_fit(sig,yr,label,obs)
end
end
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
