using CALCEPH,PyPlot
include("CGS.jl")
rc("font",family="sans-serif")
include("sim_times.jl")

function plot_orbits(dimension::Int,obs::String)
  jd1=2.4332825e6 ; nyear=10 ; sigma=30 ;jdsize=1000
  eph = Ephem("INPUTS/DE440.bsp") ; prefetch(eph)
  options = useNaifId+unitKM+unitDay # useNaifId + unitDay + unitAU
  AU = 149597870.700 #km
  jd2 = nyear*365.25 + jd1
  theta_sun=range(0, stop=2pi, length=100)
  xsun = CGS.RSUN/CGS.AU * cos.(theta_sun)
  ysun = CGS.RSUN/CGS.AU * sin.(theta_sun)
  t0 = range(jd1,stop=jd2-1,length = jdsize)
  pva_sun = zeros(6, jdsize)
  pva_venus = zeros(6, jdsize)
  pva_earth = zeros(6, jdsize)
  pva_emb = zeros(6, jdsize)
  pva_moon = zeros(6, jdsize)
  for i=1:jdsize
    pva_sun[1:6,i] = compute(eph,t0[i],0.0,10,10,options)./AU
    pva_venus[1:6,i] = compute(eph,t0[i],0.0,2,10,options)./AU
    if obs=="fromEMB"
      pva_emb[1:6,i] = compute(eph,t0[i],0.0,3,10,options) ./AU
    else
      pva_moon = compute(eph,t0[i],0.0,301,10,options)./AU
      pva_earth[1:6,i] = compute(eph,t0[i],0.0,399,10,options)./AU
    end
  end
  body,tt0,tt,sigtt=sim_obs_and_find_times(jd1,sigma,nyear,obs)
  nt1 = sum(body .== 1.0)
  nt2 = sum(body .== 2.0)
  trans_pva_venus = zeros(6, jdsize)
  trans_pva_earth = zeros(6, jdsize)
  trans_pva_emb = zeros(6, jdsize)
  trans_pva_moon = zeros(6, jdsize)
  for i=1:length(tt)
    trans_pva_venus[1:6,i] = compute(eph,tt[i],0.0,2,10,options)./AU
    if obs=="fromEMB"
      trans_pva_emb[1:6,i] = compute(eph,tt[i],0.0,3,10,options) ./AU
    else
      trans_pva_moon = compute(eph,tt[i],0.0,301,10,options)./AU
      trans_pva_earth[1:6,i] = compute(eph,tt[i],0.0,399,10,options)./AU
    end
  end
  n_obs=calc_obs_loc(trans_pva_venus[1:3],trans_pva_venus[4:6],trans_pva_earth[1:3],trans_pva_earth[4:6])
  
  if dimension==2
  fig=figure(figsize=(5,5))#,dpi=150)
  # subplot(211)
  fill(xsun,ysun,color="yellow")
  plot(xsun,ysun,color="black")
  plot(pva_venus[1,:],pva_venus[2,:],color="orange",linewidth=0.3,alpha=0.5)
  plot(pva_earth[1,:],pva_earth[2,:],color="dodgerblue",linewidth=0.3,alpha=0.5)
  scatter(trans_pva_venus[1,1:nt1],trans_pva_venus[2,1:nt1],marker=".",color="orange")
  scatter(trans_pva_earth[1,nt1+1:nt1+nt2],trans_pva_earth[2,nt1+1:nt1+nt2],marker=".",color="dodgerblue")
  plot([0,n_obs[1]*1.1],[0,n_obs[2]*1.1],"k--")
  # arrow(0.0,0.0,n_obs[1],n_obs[2],facecolor="black")
  annotate("to observer",xy=[n_obs[1];n_obs[2]], xytext=[n_obs[1]+0.05;n_obs[2]],xycoords="data") 
  xlabel("x [AU]")
  ylabel("y [AU]")
  # subplot(212)
  tight_layout()
  end

  if dimension==3
  fig=figure(figsize=(6,6))
  PyPlot.scatter3D(xsun,ysun,0,marker="o",color=:yellow)
  PyPlot.plot3D(vec(pva_venus[1,:]), vec(pva_venus[2,:]), vec(pva_venus[3,:]),alpha=0.25,color=:orange)
  PyPlot.plot3D(vec(pva_earth[1,:]), vec(pva_earth[2,:]), vec(pva_earth[3,:]),alpha=0.25,color=:dodgerblue)
  PyPlot.plot3D([0,n_obs[1]*1.2],[0,n_obs[2]*1.2],[0,n_obs[3]*1.2],linestyle="--",color=:grey)
  # PyPlot.scatter3D(vec(trans_pva_sun[1,:]), vec(trans_pva_sun[2,:]), vec(trans_pva_sun[3,:]),color=:yellow)
  PyPlot.scatter3D(vec(trans_pva_venus[1,1:nt1]),vec(trans_pva_venus[2,1:nt1]),vec(trans_pva_venus[3,1:nt1]),color=:orange,marker=".")
  PyPlot.scatter3D(vec(trans_pva_earth[1,nt1+1:nt1+nt2]),vec(trans_pva_earth[2,nt1+1:nt1+nt2]),vec(trans_pva_earth[3,nt1+1:nt1+nt2]),color=:dodgerblue,marker=".")
  PyPlot.tick_params(which="major",
      left=false,right=false,top=false,bottom=false)
  # xlim(-1,1)
  # ylim(-1,1)
  xlabel("x [AU]")
  ylabel("y [AU]")
  zlabel("z [AU]")  
  end
  #   println(n_obs) 
  # return tt 
end


function plot_res(sigma::Real,nyear::Real,sim,fitmodel,include_moon::Bool=false)
  fitfile=string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
  label="Earth"
  bestfit="best_p3"

  if String(sim)=="EMB" && isfile(string("FITS/fromEMB/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")) 
    fitfile=string("FITS/fromEMB/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
    bestfit="best_p3"
    label="EMB"
  elseif String(sim)=="EMB" && fitmodel=="p4" #if isfile(string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")) 
    fitfile=string("FITS/fromEMB",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
    bestfit="best_p4"
    label="EMB"
  elseif fitmodel=="p4" #if isfile(string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")) 
    fitfile=string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
    bestfit="best_p4"
    label="Earth"
  # else 
  #   return  println("FITS file for ",sim," with ",fitmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  f=jldopen(String(fitfile),"r")
  tt,tt0,sigtt,ttmodel=f["tt"],f["tt0"],f["sigtt"],f["ttmodel"]
  pbest_global=f[bestfit]
  nplanet,ntrans=f["nplanet"],f["ntrans"]
  # pair_ttvs=decompose_ttvs(nplanet,ntrans,f["best_p3"][1:15]) .* (24 * 60)
  p2_ttvs=decompose_ttvs(2,ntrans[1:2],f["best_p3"][1:10]) .* (24 * 60)
  p3_ttvs=decompose_ttvs(3,ntrans[1:3],f["best_p3"][1:15]) .* (24 * 60)
  n1,n2=ntrans[1],ntrans[2]
  mu1,P1,t01,ecos1,esin1=pbest_global[1:5]
  mu2,P2,t02,ecos2,esin2=pbest_global[6:10]
  time1=collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)
  time2=collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  tt1,tt2=tt[1:n1],tt[n1+1:n1+n2]
  ttsim1,ttsim2=(ttmodel[1:n1].-t01)./365.25,(ttmodel[n1+1:n1+n2].-t02)./365.25 #in years
  epoch1,epoch2=(time1.-t01)./P1,(time2.-t02)./P2
  ttv1,ttv2=(tt1.-time1).* (24 * 60),(tt2.-time2).* (24 * 60) #in minutes
  sigtt1,sigtt2=sigtt[1:n1].* (24*60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes
  # for i=1:15
  #   println("time: ",time2[i]," ttmodel-t0: ",ttsim2[i]," time-t0/Per: ",epoch2[i])
  #   # println("Venus & ",round(tt1[i],digits=5)," & ",round(ttv1[i],digits=5)," \\") 
  #   #   println(label," & ",round(tt2[i],digits=5)," & ",round(ttv2[i],digits=5)," \\") 
  # end
end