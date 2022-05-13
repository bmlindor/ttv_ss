using PyPlot,Statistics,JLD2,DelimitedFiles
rc("font",family="sans-serif")
include("decompose_ttvs.jl")
# Plot residuals to best fit TTV models for different configurations
function plot_res(sigma::Float64,nyear::Float64,sim,fitmodel,include_moon::Bool=false)
  fitfile=string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
  label="Earth"
  bestfit="best_dp"

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
  ttv1,ttv2=(tt1.-time1).* (24 * 60),(tt2.-time2).* (24 * 60) #in minutes
  sigtt1,sigtt2=sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes

  fig=figure(figsize=(8,3))
  ax1=subplot(121)
  errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")#,label="Venus")
  errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")#,label="Earth")
  plot(ttsim1,p2_ttvs[1,2,1:n1],linewidth=1,label="Venus")
  plot(ttsim2,p2_ttvs[2,1,1:n2],linewidth=1,label=label)
  ylim(-6,6)
  ylabel("TTVs [mins]")
  xlabel("Time Observed [yrs]")
  minorticks_on()
  tick_params(which="both",direction="in",top="true",right="true")
  tick_params(which="major",direction="in",top="true",right="true",length=6)
  tick_params(which="minor",direction="in",top="true",right="true",length=2)
  ax2=subplot(122)
  title("Residuals from 2-planet fit")
  plot(ttsim1,ttv1-(p2_ttvs[1,2,1:n1]))
  plot(ttsim2,ttv2-(p2_ttvs[2,1,1:n2]))
  ylim(-6,6)
  xlabel("Time [yrs]")
  ylabel("Residuals [mins]")
  minorticks_on()
  tick_params(which="both",direction="in",top="true",right="true")
  tick_params(which="major",direction="in",top="true",right="true",length=6)
  tick_params(which="minor",direction="in",top="true",right="true",length=2)
  ax1.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
      ncol=3,mode="expand",borderaxespad=0.0)
  tight_layout()

  fig=figure(figsize=(8,3))
  ax1=subplot(121)
  errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")#,label="Venus")
  errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")#,label="Earth")
  plot(ttsim1,p3_ttvs[1,3,1:n1]+p3_ttvs[1,2,1:n1],linewidth=1,label="Venus")
  plot(ttsim2,p3_ttvs[2,3,1:n2]+p3_ttvs[2,1,1:n2],linewidth=1,label=label)
  ylim(-6,6)
  ylabel("TTVs [mins]")
  xlabel("Time [yrs]")
  minorticks_on()
  tick_params(which="both",direction="in",top="true",right="true")
  tick_params(which="major",direction="in",top="true",right="true",length=6)
  tick_params(which="minor",direction="in",top="true",right="true",length=2)
  ax2=subplot(122)
 title("Residuals from 3-planet fit")
  plot(ttsim1,ttv1-(p3_ttvs[1,3,1:n1]+p3_ttvs[1,2,1:n1]))
  plot(ttsim2,ttv2-(p3_ttvs[2,3,1:n2]+p3_ttvs[2,1,1:n2]))
  ylim(-6,6)
  ylabel("Residuals [mins]")
  xlabel("Time [yrs]")
    minorticks_on()
  tick_params(which="both",direction="in",top="true",right="true")
  tick_params(which="major",direction="in",top="true",right="true",length=6)
  tick_params(which="minor",direction="in",top="true",right="true",length=2)
  ax1.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
      ncol=3,mode="expand",borderaxespad=0.0)
  tight_layout()
  if length(ntrans)==4
    p4_ttvs=decompose_ttvs(4,ntrans[1:4],f["best_p4"][1:20]) .* (24 * 60)
  fig=figure(figsize=(8,3))
  ax1=subplot(121)
    # title("4-planet fit with 10 sec noise")
    errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")#,label="Venus")
    errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")#,label="Earth")
    plot(ttsim1,p4_ttvs[1,4,1:n1]+p4_ttvs[1,3,1:n1]+p4_ttvs[1,2,1:n1],linewidth=1,label="Venus")
    plot(ttsim2,p4_ttvs[2,4,1:n2]+p4_ttvs[2,3,1:n2]+p4_ttvs[2,1,1:n2],linewidth=1,label=label)
    ylim(-6,6)
    ylabel("TTVs [mins]")
    legend(loc="upper right")
    xlabel("Time [yrs]")
    minorticks_on()
    tick_params(which="both",direction="in",top="true",right="true")
    tick_params(which="major",direction="in",top="true",right="true",length=6)
    tick_params(which="minor",direction="in",top="true",right="true",length=2)
  ax2=subplot(122)
 title("Residuals from 4-planet fit")
    plot(ttsim1,ttv1-(p4_ttvs[1,4,1:n1]+p4_ttvs[1,3,1:n1]+p4_ttvs[1,2,1:n1]))
    plot(ttsim2,ttv2-(p4_ttvs[2,4,1:n2]+p4_ttvs[2,3,1:n2]+p4_ttvs[2,1,1:n2]))
    ylim(-6,6)
    ylabel("Residuals [mins]")
    minorticks_on()
  tick_params(which="both",direction="in",top="true",right="true")
  tick_params(which="major",direction="in",top="true",right="true",length=6)
  tick_params(which="minor",direction="in",top="true",right="true",length=2)
  ax1.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
      ncol=3,mode="expand",borderaxespad=0.0)
  tight_layout()
    tight_layout()
  end
  if include_moon
    #why not loading best_p3 fit results?
  moon=moon_ttvs(ntrans[1:3],f["best_dp"]) .* (24 * 60)
  p2_ttvs=decompose_ttvs(2,ntrans[1:2],f["best_dp"][1:10]) .* (24 * 60)
  p3_ttvs=decompose_ttvs(3,ntrans[1:3],f["best_dp"][1:15]) .* (24 * 60)
  fig=figure(figsize=(8,4))
  subplot(211)
  title("3-planet + moon fit with 10 sec noise")
  errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")#,label="Venus")
  errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")#,label="Earth")
  plot(ttsim1,p3_ttvs[1,3,1:n1]+p3_ttvs[1,2,1:n1],linewidth=1,label="Venus")
  plot(ttsim2,p3_ttvs[2,3,1:n2]+p3_ttvs[2,1,1:n2]+moon,linewidth=1,label=label)
  ylabel("TTVs [mins]")
  xlabel("Time [yrs]")
  minorticks_on()
  tick_params(which="both",direction="in",top="true",right="true")
  tick_params(which="major",direction="in",top="true",right="true",length=6)
  tick_params(which="minor",direction="in",top="true",right="true",length=2)
  legend(loc="upper right")
  subplot(212)
  plot(ttsim1,ttv1-(p3_ttvs[1,3,1:n1]+p3_ttvs[1,2,1:n1]))
  plot(ttsim2,ttv2-(p3_ttvs[2,3,1:n2]+p3_ttvs[2,1,1:n2]+moon))
  ylim(-4,4)
  ylabel("Residuals [mins]")
    minorticks_on()
  tick_params(which="both",direction="in",top="true",right="true")
  tick_params(which="major",direction="in",top="true",right="true",length=6)
  tick_params(which="minor",direction="in",top="true",right="true",length=2)
  end
  tight_layout()
  # savefig("IMAGES/ttvs.eps")
  show()
end
# Plot moon signal from subtracting EMB times from Earth times
function plot_moon(sigma::Float64,nyear::Float64,fitmodel)
  EMBfit= string("FITS/fromEMB/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
  EVfit=  string("FITS/moon_fit",sigma,"s",nyear,"yrs.jld2")
  EMB=jldopen(String(EMBfit),"r")
  EV=jldopen(String(EVfit),"r")
  EMBtt,EMBtt0,EMBsigtt,EMBttmodel=EMB["tt"],EMB["tt0"],EMB["sigtt"],EMB["ttmodel"]
  EVtt,EVtt0,EVsigtt,EVttmodel=EV["tt"],EV["tt0"],EV["sigtt"],EV["ttmodel"]
  # # pair_ttvs=decompose_ttvs(nplanet,ntrans,f["best_p3"][1:15]) .* (24 * 60)
  # p2_ttvs=decompose_ttvs(2,ntrans[1:2],EMB["best_p3"][1:10]) .* (24 * 60)
  # p3_ttvs=decompose_ttvs(3,ntrans[1:3],EMB["best_p3"][1:15]) .* (24 * 60)
  EMBpbest_global=EMB["best_p3"]
  EMBnplanet,EMBntrans=EMB["nplanet"],EMB["ntrans"]
  EMBn1,EMBn2=EMBntrans[1],EMBntrans[2]
  EMBmu1,EMBP1,EMBt01,EMBecos1,EMBesin1=EMBpbest_global[1:5]
  EMBmu2,EMBP2,EMBt02,EMBecos2,EMBesin2=EMBpbest_global[6:10]
  EMBtime1=collect(EMBt01 .+ range(0,stop=EMBn1-1,length=EMBn1) .* EMBP1)
  EMBtime2=collect(EMBt02 .+ range(0,stop=EMBn2-1,length=EMBn2) .* EMBP2)
  EMBttsim1,EMBttsim2=(EMBttmodel[1:EMBn1].-EMBt01)./365.25,(EMBttmodel[EMBn1+1:EMBn1+EMBn2].-EMBt02)./365.25
  EMBtt1,EMBtt2=EMBtt[1:EMBn1],EMBtt[EMBn1+1:EMBn1+EMBn2]
  EMBttv1,EMBttv2=(EMBtt1.-EMBtime1).* (24 * 60),(EMBtt2.-EMBtime2).* (24 * 60)
  EVpbest_global=EV["best_p3"]
  EVnplanet,EVntrans=EV["nplanet"],EV["ntrans"]
  EVn1,EVn2=EVntrans[1],EVntrans[2]
  EVmu1,EVP1,EVt01,EVecos1,EVesin1=EVpbest_global[1:5]
  EVmu2,EVP2,EVt02,EVecos2,EVesin2=EVpbest_global[6:10]
  EVtime1=collect(EVt01 .+ range(0,stop=EVn1-1,length=EVn1) .* EVP1)
  EVtime2=collect(EVt02 .+ range(0,stop=EVn2-1,length=EVn2) .* EVP2)
  EVttsim1,EVttsim2=(EVttmodel[1:EVn1].-EVt01)./365.25,(EVttmodel[EVn1+1:EVn1+EVn2].-EVt02)./365.25
  EVtt1,EVtt2=EVtt[1:EVn1],EVtt[EVn1+1:EVn1+EVn2]
  EVttv1,EVttv2=(EVtt1.-EVtime1).* (24 * 60),(EVtt2.-EVtime2).* (24 * 60) #in minutes
  # sigtt1,sigtt2=sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes
  Mtt,Mtt0,Mttmodel=EMBtt.-EVtt,EMBtt0.-EVtt0,EMBttmodel.-EVttmodel
  moon=moon_ttvs(EVntrans,EV["best_dp"]) .* (24 * 60)
  Mttv=EVttv2.-EMBttv2 # transit times of actual Earth minus EMB
  fig=figure(figsize=(6,3))
  # plot(range(0,nyear,length=length(EVttv2)),EVttv2,".",linestyle="-",alpha=0.6,label="Earth signal")
  # plot(range(0,nyear,length=length(EMBttv2)),EMBttv2,".",linestyle="-",alpha=0.6,label="EMB signal")
  plot(range(0,nyear,length=length(Mttv)),Mttv,".",linestyle="-",label="Moon signal")
  plot(range(0,nyear,length=length(Mttv)),moon,label="best fit model")
  ylabel("TTVs [mins]")
  legend()
  tight_layout()
  # plot(Mtt0,Mtt,".",linestyle="-",label="Mtt")
  show()
end
# Plot timing data observed and observed - calculated (ie TTVs)
function plot_ex(sigma::Float64,nyear::Float64,sim,fitmodel,include_moon::Bool=false)
  fitfile=string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
  label="Earth"
  bestfit="best_dp"
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
  ttv1,ttv2=(tt1.-time1).* (24 * 60),(tt2.-time2).* (24 * 60) #in minutes
  sigtt1,sigtt2=sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes
  fig=plt.figure(figsize=(8,6))
  ax1=subplot(211)
  ax1.plot(ttsim1,tt1.-2430000,"o",color="grey")
  ax1.set_ylabel("Observed mid-transit time -2430000")
  ax1.set_xlabel("Time [yrs]")
  ax2=subplot(212)
  errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")
  plot(ttsim1,ttv1)
  ax2.set_ylabel("Observed - Calculated")
  ax1.set_xlabel("Time [yrs]")
    tight_layout()
  show()
end
# Plot observed TTVs vs model fit, with contributions
function plot_ttvs(sigma::Float64,nyear::Float64,sim,fitmodel,include_moon::Bool=false)
  fitfile=string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
  nplanet=3
  label="Earth contribution"
  labely="Earth TTVs [mins]"
  if String(sim)=="EMB" && isfile(string("FITS/fromEMB/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")) 
    fitfile=string("FITS/fromEMB/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
    bestfit="best_p3"
    label="EMB contribution"
    labely="EMB TTVs [mins]"
  end
  if fitmodel=="p4"
    bestfit="best_p4"
    nplanet=4
  elseif fitmodel=="p5"
    bestfit="best_p5"
    nplanet=5
  elseif fitmodel=="moon" #if isfile(string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")) 
    bestfit="best_p3"
  # else 
  #   return  println("FITS file for ",sim," with ",fitmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  f=jldopen(String(fitfile),"r")
  tt,tt0,sigtt,ttmodel=f["tt"],f["tt0"],f["sigtt"],f["ttmodel"]
  pbest_global=f[bestfit]
  ntrans=f["ntrans"]
  pair_ttvs=decompose_ttvs(nplanet,ntrans[1:nplanet],pbest_global) .* (24 * 60)
  # println(bestfit)
  n1,n2=ntrans[1],ntrans[2]
  mu1,P1,t01,ecos1,esin1=pbest_global[1:5]
  mu2,P2,t02,ecos2,esin2=pbest_global[6:10]
  # mu3,P3,t03,ecos3,esin3=pbest_global[11:15]
  time1=collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)
  time2=collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  tt1,tt2=tt[1:n1],tt[n1+1:n1+n2]
  ttsim1,ttsim2=(ttmodel[1:n1].-t01)./365.25,(ttmodel[n1+1:n1+n2].-t02)./365.25 #in years
  ttv1,ttv2=(tt1.-time1).* (24 * 60),(tt2.-time2).* (24 * 60) #in minutes
  sigtt1,sigtt2=sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes

  fig=figure(figsize=(7,5))
  ax1=subplot(211)
  # plot(ttsim1,pair_ttvs[1,3,1:n1]+pair_ttvs[1,2,1:n1],linewidth=2,color="grey")#,label="Total variations")
  # plot(ttsim1,pair_ttvs[1,3,1:n1],color="firebrick",label="Jupiter contribution",alpha=0.9)
  plot(ttsim1,pair_ttvs[1,2,1:n1],color="green",label=label,alpha=0.9)
  # plot(ttsim1,ttv1,linewidth=1.25,color="grey",label="Obs variations")
  errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")
  ylabel("Venus TTVs [mins]")
  xlabel("Time [yrs]")
  ax1.minorticks_on()
  ax1.tick_params(which="both",direction="in",top="true",right="true")
  # ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
  # ax1.tick_params(which="minor",direction="in",top="true",right="true",length=2)
  if include_moon
    moon=moon_ttvs(ntrans,pbest_global) .* (24 * 60)
    ax2=subplot(212,sharex=ax1)
    plot(ttsim2,pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2]+moon,linewidth=2,color="grey")#,label="Total variations")
    # plot(ttsim2,pair_ttvs[2,3,1:n2],color="firebrick",label="Jupiter contribution",alpha=0.9)
    plot(ttsim2,pair_ttvs[2,1,1:n2],color="orange",label="Venus contribution",alpha=0.9)
    plot(ttsim2,moon,linestyle="-.",color="purple",label="Moon contribution",alpha=0.9)
    ylabel(labely)
    errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")
    xlabel("Time [yrs]")
    ax2.minorticks_on()
    ax2.tick_params(which="both",direction="in",top="true",right="true")
  else
  ax2=subplot(212,sharex=ax1)
  # plot(ttsim2,pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2],linewidth=2,color="grey")#,label="Total variations")
  plot(ttsim2,pair_ttvs[2,1,1:n2],color="orange",label="Venus contribution",alpha=0.9)
  ylabel(labely)
  errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")
  xlabel("Time [yrs]")
  ax2.minorticks_on()
  ax2.tick_params(which="both",direction="in",top="true",right="true")
  end
  if nplanet==4
    ax1.plot(ttsim1,pair_ttvs[1,3,1:n1]+pair_ttvs[1,2,1:n1]+pair_ttvs[1,4,1:n1],linewidth=2,color="grey")
    ax1.plot(ttsim1,pair_ttvs[1,3,1:n1],linestyle="--",color="darkcyan",label="'Mars' contribution",alpha=0.9)
    ax1.plot(ttsim1,pair_ttvs[1,4,1:n1],color="firebrick",label="Jupiter contribution")
    ax2.plot(ttsim2,pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2]+pair_ttvs[2,4,1:n2],linewidth=2,color="grey")
    ax2.plot(ttsim2,pair_ttvs[2,3,1:n2],linestyle="--",color="darkcyan",label="'Mars' contribution",alpha=0.9)
    ax2.plot(ttsim2,pair_ttvs[2,4,1:n2],color="firebrick",label="Jupiter contribution")
  elseif nplanet==5
    ax1.plot(ttsim1,pair_ttvs[1,3,1:n1]+pair_ttvs[1,2,1:n1]+pair_ttvs[1,4,1:n1]+pair_ttvs[1,5,1:n1],linewidth=2,color="grey")
    ax1.plot(ttsim1,pair_ttvs[1,3,1:n1],linestyle="--",color="darkcyan",label="'Mars' contribution",alpha=0.9)
    ax1.plot(ttsim1,pair_ttvs[1,4,1:n1],color="firebrick",label="Jupiter contribution")
    ax1.plot(ttsim1,pair_ttvs[1,5,1:n1],linestyle="-.",color="tan",label="'Saturn' contribution")
    ax2.plot(ttsim2,pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2]+pair_ttvs[2,4,1:n2]+pair_ttvs[2,5,1:n2],linewidth=2,color="grey")
    ax2.plot(ttsim2,pair_ttvs[2,3,1:n2],linestyle="--",color="darkcyan",label="'Mars' contribution",alpha=0.9)
    ax2.plot(ttsim2,pair_ttvs[2,4,1:n2],color="firebrick",label="Jupiter contribution")
    ax2.plot(ttsim2,pair_ttvs[2,5,1:n2],linestyle="-.",color="tan",label="'Saturn' contribution")
  else 
    ax1.plot(ttsim1,pair_ttvs[1,3,1:n1]+pair_ttvs[1,2,1:n1],linewidth=2,color="grey")
    ax1.plot(ttsim1,pair_ttvs[1,3,1:n1],color="firebrick",label="Jupiter contribution",alpha=0.9)
    ax2.plot(ttsim2,pair_ttvs[2,3,1:n2],color="firebrick",label="Jupiter contribution",alpha=0.9)
    ax2.plot(ttsim2,pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2],linewidth=2,color="grey")
  end
  ax1.set_ylim(-6,6)
  ax2.set_ylim(-6,6)
  ax1.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
      ncol=2,mode="expand",borderaxespad=0.0)
  ax2.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
      ncol=2,mode="expand",borderaxespad=0.0)
  tight_layout()
  title=string("IMAGES/ttvs/",sim,fitmodel,"ttvs-",sigma,"secs",nyear,"yrs.png")
  # savefig(title)
  show()
end
