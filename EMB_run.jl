function show_args(args)
@show args
end
using Profile
include("sim_times.jl")
include("fit_planet3.jl")
# include("fit_planet4.jl")
# include("fit_moon.jl")
include("MCMC.jl")
show_args(ARGS)
runtype,sigma,nyear,label=ARGS[1],parse(Float64,ARGS[2]),parse(Float64,ARGS[3]),ARGS[4]
# Initialize variables and fixed values
jd1=2.4332825e6
np3,nphase,ndp=200,36,180 #wide: 100,36,180 
p3in,p3out,dpin,dpout=zeros(4)
p4in,p4out,np4= 1.5*365.25,5*365.25,100
nwalkers,nsteps=50,1000
# Change search grids to accomodate for wider distributions when time spans are shorter
if nyear>36
	p3in,p3out=11.4*365.25,12.1*365.25
	dpin,dpout=2.25,2.37
	nsteps=10000
elseif nyear>=30 && nyear<=36
	p3in,p3out=10.8*365.25,12.3*365.25
	dpin,dpout=2.2,2.38
	nsteps=50000
elseif nyear>=20 && nyear<=29
	p3in,p3out=10.6*365.25,14.2*365.25 
	dpin,dpout=2.1,2.52
	nsteps=50000
elseif nyear>=12 && nyear<=19
	p3in,p3out=10.6*365.25,14.2*365.25 
	dpin,dpout=2.1,2.52
	nsteps=50000
elseif nyear<12
	p3in,p3out=10*365.25,15*365.25 
	nsteps=100000
end
# From EMB simulations, use datafile=string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt") and EMB==true.
datafile=string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt")
p3file=string("FITS/fromEMB/p3_fit",sigma,"s",nyear,"yrs.jld2")
# Simulate orbits, create datafile, and perform 3-planet grid search
function simfit()
	sim_times(jd1,nyear,true,sigma,true)
	@time fit_planet3(datafile,jd1,sigma,nyear,p3in,p3out,np3,nphase,true,true)   
end
# Perform 3-planet grid search, given width flag
function grid_run(p3in,p3out,np3,nphase)
	@time fit_planet3(datafile,jd1,sigma,nyear,p3in,p3out,np3,nphase,true,true) 
end
# Run 3-planet markov chain
function p3_mcmc()
	foutput=string("MCMC/fromEMB/p3_mcmc",sigma,"s",nyear,"yrs.jld2")
	p=jldopen(String(p3file),"r")
	@time MCMC(foutput,p["best_p3"],p["lprob_best_p3"],nsteps,nwalkers,3,p["ntrans"][1:3],p["tt0"],p["tt"],p["sigtt"],true,true)	
end 
# Run lunar markov chain
# function moon_mcmc()
# 	foutput=string("MCMC/moon_mcmc",sigma,"s",nyear,"yrs.jld2")
# 	m=jldopen(String(p3file),"r")
# 	@time MCMC(foutput,m["best_dp"],m["lprob_best_dp"],nsteps,nwalkers,3,m["ntrans"][1:3],m["tt0"],m["tt"],m["sigtt"],true,false)	  
# end
# Run 4-planet markov chain
function p4_mcmc()
	foutput=string("MCMC/fromEMB/p4_mcmc",sigma,"s",nyear,"yrs.jld2")
	p4file=string("FITS/fromEMB/p4_fit",sigma,"s",nyear,"yrs.jld2")
	p=jldopen(String(p4file),"r")
	@time MCMC(foutput,p["best_p4"],p["lprob_best_p4"],nsteps,nwalkers,4,p["ntrans"][1:4],p["tt0"],p["tt"],p["sigtt"],true,true)
end
if runtype=="sim" && label=="ppp"
	simfit()
# elseif runtype=="mcmc" && label=="ppmp"
# 	moon_mcmc()
elseif runtype=="mcmc" && label=="ppp"
	p3_mcmc()
elseif runtype=="mcmc" && label=="pppp"
	p4_mcmc()
elseif runtype=="grid" && label=="ppp"
 	grid_run(p3in,p3out,np3,nphase)
elseif runtype=="grid" && label=="pppp" 
	@time fit_planet4(jd1,sigma,nyear,p4in,p4out,np4,nphase)
elseif runtype=="full" && label=="ppp" 
 	grid_run(p3in,p3out,np3,nphase)
 	p3_mcmc()
elseif runtype=="full" && label=="pppp"
	@time fit_planet4(jd1,sigma,nyear,p4in,p4out,np4,nphase)
	p4_mcmc()
# elseif runtype=="wide" && label=="ppmp"
# 	grid_run(3*365.25,15*365.25,200,36,0.0,2pi,180,true)
else
	println("No routine available with that runtype and/or label.")
end