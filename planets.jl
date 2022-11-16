# Julia v1.3
using CALCEPH,PyPlot,Statistics,JLD2,DelimitedFiles,Random,LinearAlgebra
rc("font",family="sans-serif")
rc("lines",linewidth=2)
include("regress.jl")
include("CGS.jl")
# Load JPL ephemerides from data and set units
eph = Ephem("INPUTS/DE440.bsp") ; prefetch(eph)
options = useNaifId+unitKM+unitDay # useNaifId + unitDay + unitAU
AU = 149597870.700 #km
Random.seed!(42)
  jd1=2.4332825e6
  sigma=30
  nyear=10
  function calc_obs_loc(pos1,vel1,pos2,vel2)
  h1 = cross(pos1,vel1)
  h2 = cross(pos2,vel2)
  n_obs = cross(h2,h1)
  n_obs /= norm(n_obs) #from one direction when both transit
end
function moon_times(jd1::Float64,sigma::Real,nyear::Real)
  # nyear = (jd2 - jd1)/365.25 
  jd2 = nyear*365.25 + jd1
  jdsize = 1000
  # dt = (jd2 - jd1)/jdsize
  @assert (jd1 >= 2287184.5) #2414105.0
  @assert (jd2 <= 2688976.5) #2488985.0
  t0 = range(jd1,stop=jd2-1,length = jdsize)

  # Compute ephemerides of Sun, Venus and Earth (or EMB):
  pva_sun = zeros(9,jdsize)
  pva_venus = zeros(9,jdsize)
  pva_earth = zeros(9,jdsize)
  #pva_emb = zeros(6, jdsize)
  pva_moon = zeros(9, jdsize)

  for i=1:jdsize
    pva_sun[1:9,i] = compute(eph,t0[i],0.0,10,10,options,2)./AU
    pva_venus[1:9,i] = compute(eph,t0[i],0.0,2,10,options,2)./AU
    pva_earth[1:9,i] = compute(eph,t0[i],0.0,399,10,options,2)./AU
      # pva_emb = compute(eph,t0[i],0.5,3,10,options,2)
    pva_moon[1:9,i] = compute(eph,t0[i],0.5,301,10,options,2)
      # println("Earth - EMB: ",norm(pva_earth[1:3,i] .- pva_emb[1:3]))
      # println("Earth - Moon: ",norm(pva_earth[1:3,i] .- pva_moon[1:3]))
      # println("Moon - EMB: ",norm(pva_moon[1:3] .- pva_emb[1:3]))
      # println("Ratio: ",norm(pva_earth[1:3,i] .- pva_emb[1:3])/norm(pva_moon[1:3] .- pva_emb[1:3]))
    end

  # Find observer location required to see transits of Venus and Earth:
  n_obs=calc_obs_loc(pva_venus[1:3],pva_venus[4:6],pva_earth[1:3],pva_earth[4:6])

  # Find actual transit times:
  P_venus = 225.0
  P_earth = 365.0
  P_err = 1.0
  tt1 = transit_times(2,eph,t0,P_venus,P_err,n_obs,10)
  nt1=length(tt1)
  tt2 = transit_times(399,eph,t0,P_earth,P_err,n_obs,10)
  nt2=length(tt2)
  tt3 = transit_times(301,eph,t0,365.25,P_err,n_obs,10)
  nt3=length(tt3)

  sigtt1=fixed_noise(tt1,sigma)
  sigtt2=fixed_noise(tt2,sigma)
  sigtt3=fixed_noise(tt3,sigma)

  x1,t01,per1 = linear_fit(tt1,P_venus,sigtt1)
  x2,t02,per2 = linear_fit(tt2,P_earth,sigtt2)
  x3,t03,per3 = linear_fit(tt3,365.25,sigtt2)
  println("moon_period=",per3)
  # println("coefficients: ",t02," , ",per2)
  tt=[tt1;tt2;tt3]

  # # Best-fit linear transit times:
  t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) 
  t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
  t3  = collect(t03 .+ per3 .* range(0,stop=nt3-1,length=nt3))
  # # tt0 = [t1;t2;t3]
  # # sigtt=[sigtt1;sigtt2;sigtt3]

  # # body = zeros((nt1+nt2))
  # # body[1:nt1] .= 1.0
  # # body[nt1+1:nt1+nt2] .= 2.0
  subplot(211)
  plot((t2.-t02)./per2,(tt2.-t2).*(24*60))
  xlabel("Time [years]")
  ylabel("TTV [min]")
  subplot(212)
  plot((t3.-t03)./per2,(tt3.-t3).*(24*60))
  xlabel("Time")
  ylabel("TTV [min]")
  toff=((tt2.-t2).*(24*60)).-((tt3.-t3).*(24*60))
  a_s=toff/2pi
  # return tt1,tt2,tt3
  # f=  jldopen(String(FITS/p3moon_fit",sigma,"s",nyear,"yrs.jld2"),"r")
  # tt,tt0,sigtt,ttmodel=f["tt"],f["tt0"],f["sigtt"],f["ttmodel"]
  # pbest_global=f["best_dp"]
  # nplanet,ntrans=f["nplanet"],f["ntrans"]
  # n1,n2=ntrans[1],ntrans[2]
  # mu1,P1,t01,ecos1,esin1=pbest_global[1:5]
  # mu2,P2,t02,ecos2,esin2=pbest_global[6:10]
  # time1=collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)  #tcalc
  # time2=collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  # tt1,tt2=tt[1:n1],tt[n1+1:n1+n2]       #tobs
  # ttmodel1,ttmodel2 = ttmodel[1:n1],ttmodel[n1+1:n1+n2]
  # ttsim1,ttsim2=(time1.-t01)./365.25,(time2.-t02)./365.25 #in years
  # ttvmodel1,ttvmodel2=(ttmodel1.-time1).*(24*60),(ttmodel2.-time2).*(24*60)
  # ttv1,ttv2=(tt1.-time1).* (24*60),(tt2.-time2).* (24 * 60) #in minutes
  # sigtt1,sigtt2=sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes
  # pair_ttvs=decompose_ttvs(f["nplanet"],ntrans[1:f["nplanet"]],pbest_global) .* (24 * 60)

  # moon=moon_ttvs(ntrans,pbest_global) .* (24 * 60)
  # subplot(211)
  # plot(ttsim2,pair_ttvs[2,1,1:n2],color="salmon",label="Venus")
  # plot(ttsim2,moon,linestyle="-.",color="purple",label="Moon")
  # errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",capsize=3,ms=5)#,label="Earth")
  # subplot(212)
  # plot(ttsim3,moon,linestyle="-.",color="purple",label="Moon")
  # errorbar(ttsim3,ttv3,sigtt3,fmt=".",color="black",capsize=3,ms=5)#,label="Earth")
  # # text(0,-5.5,label,fontsize="xx-large")
  # xlabel("Time [years]",fontsize=20)
  # ylabel("TTV [min]",fontsize=20)
  # ylim(-7,7)
  # minorticks_on()
  # tick_params(which="both",direction="in",top=true,right=true)
end
moon_times(jd1,sigma,nyear)
