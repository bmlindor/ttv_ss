using Profile,TTVFaster,DelimitedFiles,JLD2,DataFrames,CSV,LsqFit,Statistics
include("regress.jl")
include("sim_times.jl")
include("fit_planet2.jl")
include("fit_planet3.jl")
include("fit_planet4.jl")
include("fit_planet5.jl")
include("fit_moon.jl")
include("MCMC.jl")
label,runtype,obs,opt=ARGS[1],ARGS[2],ARGS[3],parse(Int64,ARGS[4])
sigma,nyear=parse(Int64,ARGS[5]),parse(Int64,ARGS[6])
function show_args(args)
@show args
end
function parse_model(label::String)
	nplanet = 0
	nmoon = 0
	for obj in 1:length(label)
	  if label[obj]=='p'
	    nplanet += 1
	  end
	  if label[obj]=='m'
	  	nmoon += 1
	  end
	end
	println("Modeling ",nplanet," planets and ", nmoon," moons.")
	return nplanet,nmoon
end
nplanet,nmoon=parse_model(label)
show_args(ARGS)
# Initialize variables and period ranges
jd1=2.4332825e6
tref=2430000; tol=1e-5
nphase=36 #wide: 100,36,180
if opt < 1000 
nper=opt
nsteps=10000
else
nsteps=opt
nper=100
end
np3,np4,ndp,np5=[nper, nper, nper, nper]

# Run markov chains
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
	println(p," loaded.")
	MCMC(foutput,p[string("best_p",nplanet)],p[string("lprob_best_p",nplanet)],nsteps,nwalkers,nplanet,p["ntrans"][1:nplanet],p["tt0"],p["tt"],p["sigtt"],true,true)
end 
function moon_mcmc(sigma,nyear,nplanet,nsteps,label::String)
	if label=="Hppmpp" # load p4 search after moon fit
	fitfile=string("FITS/p3moonp4_fit",sigma,"s",nyear,"yrs.jld2")
	foutput=string("MCMC/p3moonp4_mcmc",sigma,"s",nyear,"yrs.jld2")
	m=jldopen(String(fitfile),"r")
	MCMC(foutput,m["best_p4"],m["lprob_best_p4"],nsteps,nwalkers,4,m["ntrans"][1:4],m["tt0"],m["tt"],m["sigtt"],true,false)
	else
	fitfile=string("FITS/p",nplanet,"moon_fit",sigma,"s",nyear,"yrs.jld2")
	foutput=string("MCMC/p",nplanet,"moon_mcmc",sigma,"s",nyear,"yrs.jld2")
	m=jldopen(String(fitfile),"r")
	println(m," loaded.")
	MCMC(foutput,m["best_dp"],m["lprob_best_dp"],nsteps,nwalkers,nplanet,m["ntrans"][1:nplanet],m["tt0"],m["tt"],m["sigtt"],true,false)	
	end
end
# Perform grid searches
# Change search grids to accomodate for wider distributions when time spans are shorter

function run_grid(sigma,nyear,label::String,obs::String)
 if nyear>=26
 	p3in,p3out=11*365.25,13*365.25
	p4in,p4out=1.8*365.25,1.95*365.25
	p5in,p5out=28*365.25,30*365.25
 	dpin,dpout=2.25,2.37
 elseif 22<=nyear<26
 	p3in,p3out=11*365.25,13*365.25
	p4in,p4out=1.8*365.25,2*365.25
	p5in,p5out=28*365.25,32*365.25
 	dpin,dpout=2.2,2.38
 elseif 18<=nyear<22
 	p3in,p3out=11*365.25,13*365.25
	p4in,p4out=1.8*365.25,2.2*365.25
	p5in,p5out=26*365.25,34*365.25
 	dpin,dpout=2.2,2.38
 elseif 12<=nyear<18
 	p3in,p3out=11*365.25,13*365.25 
	p4in,p4out=1.8*365.25,3.6*365.25
 	dpin,dpout=2.1,2.52
 	nsteps=50000
 elseif nyear<12
 	p3in,p3out=11*365.25,13*365.25 
	p4in,p4out=1.6*365.25,3.8*365.25
 end
	if label=="Hpp" 
	sim_times(jd1,sigma,nyear,obs)
 	fit_planet2(jd1,sigma,nyear,tref,tol,[obs],true)
	elseif label=="Hppp"
	fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,[obs,"p3"],true)
	elseif label=="Hpppp"
  	fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,[obs,"p4"],true)
	elseif label=="Hppppp"
	fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,[obs,"p5"],true)
	end
	if label=="Hppm"
 	#fit_planet2(jd1,sigma,nyear,tref,tol,[obs],true)
	fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,["p2moon"],true,2)
	elseif label=="Hppmp"
	#fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,[obs,"p3"],true)
	fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,["p3moon"],true)
	elseif label=="Hppmpp"
  	#fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,[obs,"p4"],true)
	# fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,4,"p4moon",true)
	fit_moon(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,["p3moonp4"],true)
	elseif label=="Hppmppp"
  	fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,[obs,"p5"],true)
	fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,["p5moon"],true,5)
	end
end

function wide_run(sigma,nyear,label,obs)
	nphase=36 #wide: 100,36,180 
	p3in,p3out,np3=5*365.25,22*365.25,nper
	dpin,dpout,ndp=0.0,2*pi,nper
	p4in,p4out,np4=1.5*365.25,5*365.25,nper
	p5in,p5out,np5=22*365.25,36*365.25,nper
	if label=="Hpp" 
	sim_times(jd1,sigma,nyear,obs)
 	best_p2,lprob_best_p2=fit_planet2(jd1,sigma,nyear,tref,tol,[obs],true)
	elseif label=="Hppp"
	fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,[obs,"widep3"],true)
	elseif label=="Hpppp"
  	fit_planet4(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,[obs,"widep4"],true)
	elseif label=="Hppppp"
  	fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,[obs,"widep5"],true)
	elseif label=="Hppmpp" 
	fit_moon(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,["widep3moonp4"],true)
	elseif label=="Hppmp"
	fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,["widep3moon"],true)

	end
end

nyears=[30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10]
#if nplanet==2 || nplanet==3
#nyears=[30,29,28,27,26,25,24,23,22,21,20,19,18]
sigmas=[10,60]#,80,90]
#end

nwalkers=75
#nsteps=50000
#for sig in sigmas
#for yr in nyears
if runtype=="mcmc" && nmoon==0
	@time planet_mcmc(sigma,nyear,nplanet,nsteps,obs)
	elseif runtype=="mcmc" && nmoon>0
	@time moon_mcmc(sigma,nyear,nplanet,nsteps,label)
end
if runtype=="grid"
	@time run_grid(sigma,nyear,label,obs)
end
if runtype=="wide"
	@time wide_run(sigma,nyear,label,obs)
end
if runtype=="widemcmc" 
	if nmoon==0
		fitfile=string("FITS/widep",nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
       	foutput=string("MCMC/widep",nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
       	p=jldopen(String(fitfile),"r")
      	println(p," loaded.")
        @time MCMC(foutput,p[string("best_p",nplanet)],p[string("lprob_best_p",nplanet)],nsteps,nwalkers,nplanet,p["ntrans"][1:nplanet],p["tt0"],p["tt"],p["sigtt"],true,true)
    elseif nmoon>0 && nplanet>3
    	fitfile=string("FITS/widep",nplanet-1,"moonp",nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
       	foutput=string("MCMC/widep",nplanet-1,"moonp",nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
       	p=jldopen(String(fitfile),"r")
      	println(p," loaded.")
        @time MCMC(foutput,p[string("best_dp")],p[string("lprob_best_dp")],nsteps,nwalkers,nplanet,p["ntrans"][1:nplanet],p["tt0"],p["tt"],p["sigtt"],true,false)
   end
end

