using Test,CALCEPH,PyPlot
include("sim_times.jl")
eph = Ephem("/Users/bethleelindor/work/washington/premap2022/jup344.bsp") ; prefetch(eph)
options = useNaifId+unitKM+unitDay # useNaifId + unitDay + unitAU
AU = 149597870.700 #km
jd1=2.4332825e6
tref=2430000; tol=1e-5
jdsize=10
sigma=30
nyear=10
obs="fromEMB"
nplanets=4
dir="test" 
dir_name=joinpath(@__DIR__,dir) ## set directory name
if !ispath(dir_name)
  mkpath(dir_name)
  return dir_name
end
function test_sim_obs_and_find_times(sigma,nyear)
 sim_times(jd1,sigma,nyear,obs)
end
# OUTPUT: .txt file with body_number (n), initial time of transit (tt0), actual transit times (tt), measurement error (sigtt)
function test_fit(sigma,nyear,nplanets,nmoon::Real=0)
  filename=string(dir,"/tt_",sigma,"s",nyear,"yrs.txt")
  jmax=5
  nphase=36 #wide: 100,36,180 
  p3in,p3out,np3=11.4*365.25,12.2*365.25,20
  p4in,p4out,np4=1.6*365.25,3*365.25,20
  p5in,p5out,np5=28*365.25,30*365.25,20
  dpin,dpout,ndp=2.29,2.35,10
  if nplanets==2
  fit_planet2(filename,jmax,jd1,sigma,nyear,tref,tol,obs,dir)
  elseif nplanets==3
  fit_planet3(filename,jmax,jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,obs,dir)
  elseif nplanets==4
  fit_planet4(filename,jmax,jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,obs,dir)
  elseif nplanets==5
  fit_planet5(filename,jmax,jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,obs,dir)
  # If planet fit exists, can do moon search
  elseif nmoon==1 && isfile(string(dir,"/p",nplanets,"_fit",sigma,"s",nyear,"yrs.jld2"))
  fit_moon(jd1,sigma,nyear,tref,tol,dpin,dpout,ndp,[string(nplanets,"moon")],nplanets)
  # elseif nmoon==1 && nplanets==4 
  # fit_moon(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,["p3moonp4"],true)
  end
end 
# OUTPUT: .jld2 file with best_p{nplanets} lprob_best_p{nplanets} ntrans nplanet tt0 tt ttmodel sigtt
function test_mcmc(nplanets,sigma,nyear)
  nwalkers=50
  nsteps=1000
  fitfile=string(dir,"/p",nplanets,"_fit",sigma,"s",nyear,"yrs.jld2")
  foutput=string(dir,"/p",nplanets,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  p=jldopen(String(fitfile),"r")
  MCMC(foutput,p[string("best_p",nplanets)],p[string("lprob_best_p",nplanets)],nsteps,nwalkers,nplanets,p["ntrans"][1:nplanets],p["tt0"],p["tt"],p["sigtt"],true,true)
  m=jldopen(foutput,"r")
  par_mcmc=m["par_mcmc"]
  iburn=m["iburn"]
  med_params=[median(par_mcmc[:,iburn:end,i]) for i=1:length(par_mcmc[1,end,:])]
  BIC(m["tt0"],nplanets,m["ntrans"],med_params,m["tt"],m["sigtt"],jmax,true)
  return
end

end
# body,tt0,tt,sigtt=test_sim_obs(jd1,sigma,nyear,obs)

# JD_venus,ff_venus,i_min_venus,pos_venus,tt_venus = find_transit(2,eph,jd1,jd1+365,n_obs,365)
# JD_earth,ff_earth,i_min_earth,pos_earth,tt_earth = find_transit(3,eph,jd1,jd1+365,n_obs,365)

# P_venus = 225.0
# P_earth = 365.0
# P_err = 1.0
# tt1 = transit_times(2,eph,t0,P_venus,P_err,n_obs,10)
# tt2 = transit_times(3,eph,t0,P_earth,P_err,n_obs,10)
# subplot(211)
# plot(JD_venus,ff_venus)
# plot(JD_earth,ff_earth)
# scatter(tt_venus,tt_venus)
# scatter(tt_earth,tt_earth)
# subplot(212)
# plot(,tt_earth)
# plot(,tt_venus)
# show()
# @testset "TTV_SS" begin
#     @test test_sim_obs_and_find_times()
#     @test test_fit()
#     @test test_mcmc()
# end