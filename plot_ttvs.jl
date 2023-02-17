using TTVFaster,PyPlot,Statistics,JLD2,DelimitedFiles
rc("font",family="sans-serif")
rc("lines",linewidth=2)
include("decompose_ttvs.jl")
include("histogram.jl")
# Plot residuals to best fit TTV models for different configurations
function plot_res(sigma::Real,nyear::Real,options,include_moon::Bool=false)
  obs=options[1]; fit_type_nplanet=options[2]; bestfit=options[3]
  if obs=="fromEMB"
    fitfile=string("FITS/fromEMB/",fit_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    fitfile2=string("FITS/fromEMB/p2_fit",sigma,"s",nyear,"yrs.jld2")
    fitfile3=string("FITS/fromEMB/p3_fit",sigma,"s",nyear,"yrs.jld2")
    label="EMB"
  elseif obs=="fromEV"
    fitfile=string("FITS/",fit_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    fitfile2=string("FITS/p2_fit",sigma,"s",nyear,"yrs.jld2")
    fitfile3=string("FITS/p3_fit",sigma,"s",nyear,"yrs.jld2")
    fitfile4=string("FITS/p4_fit",sigma,"s",nyear,"yrs.jld2")
    f4=jldopen(String(fitfile4),"r")
    label="Earth"
  #   return  println("FITS file for ",sim," with ",fitmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  f=jldopen(String(fitfile),"r")
  f2=jldopen(String(fitfile2),"r")
  f3=jldopen(String(fitfile3),"r")
  tt,tt0,sigtt,ttmodel=f["tt"],f["tt0"],f["sigtt"],f["ttmodel"]
  pbest_global=f[bestfit]
  nplanet,ntrans=f["nplanet"],f["ntrans"]
  # pair_ttvs=decompose_ttvs(nplanet,ntrans,f["best_p3"][1:15]) .* (24 * 60)
  p2_ttvs=decompose_ttvs(2,ntrans[1:2],f2["best_p2"][1:10]) .* (24 * 60)
  p3_ttvs=decompose_ttvs(3,ntrans[1:3],f3["best_p3"][1:15]) .* (24 * 60)
  n1,n2=ntrans[1],ntrans[2]
  mu1,P1,t01,ecos1,esin1=pbest_global[1:5]
  mu2,P2,t02,ecos2,esin2=pbest_global[6:10]
  time1=collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)
  time2=collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  tt1,tt2=tt[1:n1],tt[n1+1:n1+n2]
  ttmodel1,ttmodel2 = ttmodel[1:n1],ttmodel[n1+1:n1+n2]
  ttsim1,ttsim2=(time1.-t01)./365.25,(time2.-t02)./365.25 #in years
  ttvmodel1,ttvmodel2=(ttmodel1.-time1).*(24*60),(ttmodel2.-time2).*(24*60)
  ttv1,ttv2=(tt1.-time1).* (24 * 60),(tt2.-time2).* (24 * 60) #in minutes
  sigtt1,sigtt2=sigtt[1:n1].* (24*60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes

  fig=figure(figsize=(8,6),dpi=150)
    subplots_adjust(hspace=0.25,wspace=0.05)
  # subplots_adjust(hspace=0.05,wspace=0.05)
    # title("Test")
  # title="TTVFaster Models"
  # suptitle(string(title," [",nyear," yr span], ",L"$\sigma_{obs}=$",sigma," sec "))
  ax1=subplot(221)
  ax1.errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")#,label="Venus")
  ax1.plot(ttsim1,p2_ttvs[1,2,1:n1],label=L"$\mathcal{H}_{PP}$")
  ax1.plot(ttsim1,p3_ttvs[1,3,1:n1]+p3_ttvs[1,2,1:n1],label=L"$\mathcal{H}_{PPP}$")
  ax1.text(-1,6.5,"Venus",fontweight="bold",fontsize="medium")
  ylim(-8,8)
  xlabel("Time [yrs]")
  ylabel(L"$O - C$ [mins]")
  minorticks_on()
  tick_params(which="both",direction="in",top=true,right=true)
  ax2=subplot(222)
  ax2.axhline(y=0, color="grey",linewidth=1.0,linestyle="--")
  ax2.plot(ttsim1,ttv1-(p2_ttvs[1,2,1:n1]),label=L"$\mathcal{H}_{PP}$")
  ax2.plot(ttsim1,ttv1-(p3_ttvs[1,3,1:n1]+p3_ttvs[1,2,1:n1]),label=L"$\mathcal{H}_{PPP}$")
  ylim(-8,8)
  xlabel("Time [yrs]")
  ylabel("Residuals [mins]")
  minorticks_on()
  tick_params(which="both",direction="in",top=true,right=true)

  # fig=figure(figsize=(8,3))
  ax3=subplot(223)
  ax3.errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")#,label="Earth")
  ax3.plot(ttsim2,p2_ttvs[2,1,1:n2],label=L"$\mathcal{H}_{PP}$")
  ax3.plot(ttsim2,p3_ttvs[2,3,1:n2]+p3_ttvs[2,1,1:n2],label=L"$\mathcal{H}_{PPP}$")
  ax3.text(-1,6.5,label,fontweight="bold",fontsize="medium")
  ylim(-8,8)
  ylabel(L"$O - C$ [mins]")
  xlabel("Time [yrs]")
  minorticks_on()
  tick_params(which="both",direction="in",top=true,right=true)
  ax4=subplot(224)
  ax4.axhline(y=0, color="grey",linewidth=1.0,linestyle="--")
  ax4.plot(ttsim2,ttv2-(p2_ttvs[2,1,1:n2]),label=L"$\mathcal{H}_{PP}$")
  ax4.plot(ttsim2,ttv2-(p3_ttvs[2,3,1:n2]+p3_ttvs[2,1,1:n2]),label=L"$\mathcal{H}_{PPP}$")
  ylim(-8,8)
  ylabel("Residuals [mins]")
  xlabel("Time [yrs]")
  minorticks_on()
  tick_params(which="both",direction="in",top=true,right=true)

  # ax1.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
  #     ncol=3,mode="expand",borderaxespad=0.0)

  # if nplanet==4
  #   p4_ttvs=decompose_ttvs(4,ntrans[1:4],f4["best_p4"][1:20]) .* (24 * 60)
  #   ax1.plot(ttsim1,p4_ttvs[1,4,1:n1]+p4_ttvs[1,3,1:n1]+p4_ttvs[1,2,1:n1],linestyle="-.",label=L"$\mathcal{H}_{PPPP}$")
  #   ax2.plot(ttsim1,ttv1-(p4_ttvs[1,4,1:n1]+p4_ttvs[1,3,1:n1]+p4_ttvs[1,2,1:n1]),linestyle="-.",label=L"$\mathcal{H}_{PPPP}$")
  #   ax3.plot(ttsim2,p4_ttvs[2,4,1:n2]+p4_ttvs[2,3,1:n2]+p4_ttvs[2,1,1:n2],linestyle="-.",label=L"$\mathcal{H}_{PPPP}$")
  #   ax4.plot(ttsim2,ttv2-(p4_ttvs[2,4,1:n2]+p4_ttvs[2,3,1:n2]+p4_ttvs[2,1,1:n2]),linestyle="-.",label=L"$\mathcal{H}_{PPPP}$")
  #   # ax2.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
  #   #     ncol=3,mode="expand",borderaxespad=0.0)
  # end
  if include_moon
    #why not loading best_p3 fit results?
  moon3=moon_ttvs(ntrans[1:3],f["best_dp"]) .* (24 * 60)
  p3_ttvs=decompose_ttvs(3,ntrans[1:3],f["best_dp"][1:15]) .* (24 * 60)
  # fig=figure(figsize=(8,4))
  # subplot(211)
  # errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")#,label="Venus")
  # errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")#,label="Earth")
  ax1.plot(ttsim1,p3_ttvs[1,3,1:n1]+p3_ttvs[1,2,1:n1],linestyle="--",label=L"$\mathcal{H}_{PPsP}$")
  ax3.plot(ttsim2,p3_ttvs[2,3,1:n2]+p3_ttvs[2,1,1:n2]+moon3,linestyle="--",label=L"$\mathcal{H}_{PPsP}$")
  ax2.plot(ttsim1,ttv1-(p3_ttvs[1,3,1:n1]+p3_ttvs[1,2,1:n1]),linestyle="--",label=L"$\mathcal{H}_{PPsP}$")
  ax4.plot(ttsim2,ttv2-(p3_ttvs[2,3,1:n2]+p3_ttvs[2,1,1:n2]+moon3),linestyle="--",label=L"$\mathcal{H}_{PPsP}$")
  moon4=moon_ttvs(ntrans[1:4],f["best_p4"]) .* (24 * 60)
  p4_ttvs=decompose_ttvs(4,ntrans[1:4],f["best_p4"][1:20]) .* (24 * 60)
  ax1.plot(ttsim1,p4_ttvs[1,4,1:n1]+p4_ttvs[1,3,1:n1]+p4_ttvs[1,2,1:n1],linestyle="-.",label=L"$\mathcal{H}_{PPsPP}$")
  ax3.plot(ttsim2,p4_ttvs[2,4,1:n2]+p4_ttvs[2,3,1:n2]+p4_ttvs[2,1,1:n2]+moon4,linestyle="-.",label=L"$\mathcal{H}_{PPsPP}$")
  ax2.plot(ttsim1,ttv1-(p4_ttvs[1,4,1:n1]+p4_ttvs[1,3,1:n1]+p4_ttvs[1,2,1:n1]),linestyle="-.",label=L"$\mathcal{H}_{PPsPP}$")
  ax4.plot(ttsim2,ttv2-(p4_ttvs[2,4,1:n2]+p4_ttvs[2,3,1:n2]+p4_ttvs[2,1,1:n2]+moon4),linestyle="-.",label=L"$\mathcal{H}_{PPsPP}$")  
end
  ax2.legend(loc="upper right",ncol=3)
  ax4.legend(loc="upper right",ncol=3)
  tight_layout()
  # savefig("IMAGES/EVresiduals.eps")
  show()
end
# Plot moon signal from subtracting EMB times from Earth times
function plot_moon(sigma::Real,nyear::Real,fitmodel)
  EMBfit= string("FITS/fromEMB/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
  EVfit=  string("FITS/p3moon_fit",sigma,"s",nyear,"yrs.jld2")
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
  EVpbest_global=EV["best_dp"][1:15]
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
  fig=figure(figsize=(8,3))
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
function plot_ttv(sigma::Real,nyear::Real,options::Array{String},include_moon::Bool=false)
	obs=options[1]; fit_type_nplanet=options[2]; bestfit=options[3]
  if obs=="fromEMB"
    fitfile=string("FITS/fromEMB/",fit_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    label="EMB"
  elseif obs=="fromEV"
  	fitfile=string("FITS/",fit_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
  	label="Earth"
  #   return  println("FITS file for ",sim," with ",fitmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  f=jldopen(String(fitfile),"r")
  tt,tt0,sigtt,ttmodel=f["tt"],f["tt0"],f["sigtt"],f["ttmodel"]
  pbest_global=f[bestfit]
  nplanet,ntrans=f["nplanet"],f["ntrans"]
  n1,n2=ntrans[1],ntrans[2]
  mu1,P1,t01,ecos1,esin1=pbest_global[1:5]
  mu2,P2,t02,ecos2,esin2=pbest_global[6:10]
  time1=collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1) 	#tcalc
  time2=collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  tt1,tt2=tt[1:n1],tt[n1+1:n1+n2] 			#tobs
  ttmodel1,ttmodel2 = ttmodel[1:n1],ttmodel[n1+1:n1+n2]
  ttsim1,ttsim2=(time1.-t01)./365.25,(time2.-t02)./365.25 #in years
	ttvmodel1,ttvmodel2=(ttmodel1.-time1).*(24*60),(ttmodel2.-time2).*(24*60)
  ttv1,ttv2=(tt1.-time1).* (24*60),(tt2.-time2).* (24 * 60) #in minutes
  sigtt1,sigtt2=sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes
    
  fig=plt.figure(figsize=(8,6))
  subplots_adjust(hspace=0.0)
  ax1=subplot(211) 
  title("2-planet fit of Venus and Earth",fontsize="x-large")
	ax1.minorticks_on()
  #setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels
	plot(ttsim1,ttvmodel1,label="2-planet model")
  errorbar(ttsim1,ttv1,sigtt1,fmt="v",color="black",capsize=3,ms=3)
  tick_params(which="both",direction="in",top=true,right=true)
  #plot(ttsim1,p2_ttvs[1,2,1:n1],color="salmon")
	text(0,5,"Venus",fontsize="xx-large")
  ylabel("TTV [min]",fontsize="x-large")
	ylim(-7,7)
  legend(loc="lower left")


  ax2=subplot(212,sharex=ax1) # Create the 1st axis of a 3x1 array of axes
  plot(ttsim2,ttvmodel2,label="2-planet model")
  errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",capsize=3,ms=5)#,label="Earth")
  #plot(ttsim2,p2_ttvs[2,1,1:n2],color="darkcyan")
	text(0,-5.5,label,fontsize="xx-large")
  ylabel("TTV [min]",fontsize="x-large")
  xlabel("Time [years]",fontsize="x-large")
  ylim(-7,7)
	ax2.minorticks_on()
  # yticks(fontsize=12)
  tick_params(which="both",direction="in",top=true,right=true)
  sim_obs_label= string(L"$\sigma_{obs}=$",sigma," sec")
  ax1.text(maximum(ttsim1)-4,4,sim_obs_label)
  # tick_params(which="major",direction="in",top=true,right=true,length=6)
  # tick_params(which="minor",direction="in",top=true,right=true,length=2)
  # yticks(0.1:0.2:0.9) # Set the y-tick range and step size, 0.1 to 0.9 in increments of 0.2
  # ylim(0.0,1.0) # Set the y-limits from 0.0 to 1.0

  # ax1=subplot(211)
  # # ax1.plot(ttsim1,tt1.-2430000,"o",color="grey")
  # errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")
  # plot(ttsim1,ttv1)
  # ax1.set_ylabel("Observed mid-transit time -2430000")
  # ax1.set_xlabel("Time [yrs]")
  # ax2=subplot(212)
  # errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")
  # plot(ttsim1,ttv1)
  # ax2.set_ylabel("Observed - Calculated")
  # ax1.set_xlabel("Time [yrs]")
  #show()
end
# Plot observed TTVs vs model fit, with contributions
function plot_contrib(sigma::Real,nyear::Real,options::Array{String},include_moon::Bool=false)
  obs=options[1]; fit_type_nplanet=options[2]; bestfit=options[3]
  if obs=="fromEMB"
    fitfile=string("FITS/fromEMB/",fit_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    label="EMB"
  elseif obs=="fromEV"
    fitfile=string("FITS/",fit_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    label="Earth"
  #   return  println("FITS file for ",sim," with ",fitmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  f=jldopen(String(fitfile),"r")
  tt,tt0,sigtt,ttmodel=f["tt"],f["tt0"],f["sigtt"],f["ttmodel"]
  pbest_global=f[bestfit]
  ntrans=f["ntrans"]
  nplanet=f["nplanet"]
  pair_ttvs=decompose_ttvs(nplanet,ntrans[1:nplanet],pbest_global) .* (24 * 60)
  # println(bestfit)
  n1,n2=ntrans[1],ntrans[2]
  mu1,P1,t01,ecos1,esin1=pbest_global[1:5]
  mu2,P2,t02,ecos2,esin2=pbest_global[6:10]
  # mu3,P3,t03,ecos3,esin3=pbest_global[11:15]
  time1=collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)
  time2=collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  tt1,tt2=tt[1:n1],tt[n1+1:n1+n2]
  ttmodel1,ttmodel2 = ttmodel[1:n1],ttmodel[n1+1:n1+n2]
  ttsim1,ttsim2=(time1.-t01)./365.25,(time2.-t02)./365.25 #in years
  ttvmodel1,ttvmodel2=(ttmodel1.-time1).*(24*60),(ttmodel2.-time2).*(24*60)
  ttv1,ttv2=(tt1.-time1).* (24 * 60),(tt2.-time2).* (24 * 60) #in minutes
  sigtt1,sigtt2=sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) 
  total1=0;total2=0

  fig=figure(figsize=(5,5),dpi=150)
  # title="Planet Contributions"
  # suptitle(string(title," [",nyear," yr span], ",L"$\sigma_{obs}=$",sigma," sec "))
  subplots_adjust(hspace=0.25)
  ax1=subplot(211)
  plot(ttsim1,pair_ttvs[1,2,1:n1],color="forestgreen",label=label)
  errorbar(ttsim1,ttv1,sigtt1,fmt="v",color="black",mfc="white",capsize=3,ms=3)
  ax1.text(-1,5,"Venus",fontweight="bold",fontsize="medium")
  # text(0,5,"Venus",fontsize="xx-large")
  ylabel("TTV [min]",fontsize="x-large")
  xlabel("Time [years]",fontsize="x-large")
  ylim(-8,8)
  ax1.minorticks_on()
  tick_params(which="both",direction="in",top=true,right=true)
  ax2=subplot(212,sharex=ax1)
  plot(ttsim2,pair_ttvs[2,1,1:n2],color="salmon",label="Venus")
  errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mfc="white",capsize=3,ms=5)#,label="Earth")
  ax2.text(-1,6,label,fontweight="bold",fontsize="medium")
  # text(0,-5.5,label,fontsize="xx-large")
  xlabel("Time [years]",fontsize="x-large")
  ylabel("TTV [min]",fontsize="x-large")
  ylim(-8,8)
  ax2.minorticks_on()
  tick_params(which="both",direction="in",top=true,right=true)
  if f["nplanet"]==4
    total1=pair_ttvs[1,3,1:n1]+pair_ttvs[1,2,1:n1]+pair_ttvs[1,4,1:n1]
    total2=pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2]+pair_ttvs[2,4,1:n2]
    ax1.plot(ttsim1,pair_ttvs[1,3,1:n1],linestyle="--",color="orange",label="Mars",alpha=0.9)
    ax1.plot(ttsim1,pair_ttvs[1,4,1:n1],color="firebrick",label="Jupiter")
    ax2.plot(ttsim2,pair_ttvs[2,3,1:n2],linestyle="--",color="orange",label="Mars",alpha=0.9)
    ax2.plot(ttsim2,pair_ttvs[2,4,1:n2],color="firebrick",label="Jupiter")
    ax1.set_title(L"Contributions to $\mathcal{H}_{PPPP}$ Fit",fontsize="x-large")
  elseif f["nplanet"]==5
    total1=pair_ttvs[1,3,1:n1]+pair_ttvs[1,2,1:n1]+pair_ttvs[1,4,1:n1]+pair_ttvs[1,5,1:n1]
    total2=pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2]+pair_ttvs[2,4,1:n2]+pair_ttvs[2,5,1:n2]
    ax1.plot(ttsim1,pair_ttvs[1,3,1:n1],linestyle="--",color="orange",label="'Mars'",alpha=0.9)
    ax1.plot(ttsim1,pair_ttvs[1,4,1:n1],color="firebrick",label="Jupiter")
    ax1.plot(ttsim1,pair_ttvs[1,5,1:n1],linestyle="--",color="tan",label="'Saturn'")
    ax2.plot(ttsim2,pair_ttvs[2,3,1:n2],linestyle="--",color="orange",label="'Mars'",alpha=0.9)
    ax2.plot(ttsim2,pair_ttvs[2,4,1:n2],color="firebrick",label="Jupiter")
    ax2.plot(ttsim2,pair_ttvs[2,5,1:n2],linestyle="--",color="tan",label="'Saturn'")
    ax1.set_title(L"Contributions to $\mathcal{H}_{PPPP}$ Fit",fontsize="x-large")
  elseif f["nplanet"]==3
    total1=pair_ttvs[1,3,1:n1]+pair_ttvs[1,2,1:n1]
    total2=pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2]
    ax1.plot(ttsim1,pair_ttvs[1,3,1:n1],color="firebrick",label="Jupiter")
    ax2.plot(ttsim2,pair_ttvs[2,3,1:n2],color="firebrick",label="Jupiter")
    ax1.set_title(L"Contributions to $\mathcal{H}_{PPP}$ Fit",fontsize="x-large")
  end
  ax1.plot(ttsim1,total1,linewidth=1.5,color="grey")
  if include_moon 
    moon=moon_ttvs(ntrans,pbest_global) .* (24 * 60)
    ax2.plot(ttsim2,total2+moon,linewidth=1.5,color="grey")
    ax2.plot(ttsim2,moon,linestyle="-.",color="purple",label="Moon")
    # text(0,-5.5,label,fontsize="xx-large")
    # xlabel("Time [years]",fontsize=20)
    # ylabel("TTV [min]",fontsize=20)
    # ylim(-7,7)
    # ax2.minorticks_on()
    # tick_params(which="both",direction="in",top=true,right=true)
  else
    ax2.plot(ttsim2,total2,linewidth=1.5,color="grey")
  end
  sim_obs_label= string(L"$\sigma_{obs}=$",sigma," sec")
  # ax1.set_title(L"Contributions to $\mathcal{H}_{PPPP}$ Fit",fontsize="x-large")
  ax1.text(maximum(ttsim1)-5,5,sim_obs_label)
  # ax2.text(maximum(ttsim1)-5,5,sim_obs_label)
  ax1.legend(loc="lower left",fontsize="large",ncol=nplanet)
  # ax1.legend(loc="lower left",fontsize="large",title="Contributions",title_fontsize="large",ncol=4)
  ax2.legend(loc="lower left",ncol=nplanet,fontsize="medium")
  tight_layout()
  # legend(loc="upper right")
  # title=string("IMAGES/ttvs/",sim,fitmodel,"ttvs-",sigma,"secs",nyear,"yrs.png")
  # savefig(title)
  # show()
end
  #how does scatter in residuals compare w/ uncertainty?
function plot_scatter(sigma::Real,nyear::Real,options::Array{String},include_moon::Bool=false)
  obs=options[1]; fit_type_nplanet=options[2]; bestfit=options[3]
  if obs=="fromEMB"
    fitfile=string("FITS/fromEMB/",fit_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    label="EMB"
  elseif obs=="fromEV"
    fitfile=string("FITS/",fit_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    label="Earth"
  #   return  println("FITS file for ",sim," with ",fitmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  f=jldopen(String(fitfile),"r")
  tt,tt0,sigtt,ttmodel=f["tt"],f["tt0"],f["sigtt"],f["ttmodel"]
  pbest_global=f[bestfit]
  ntrans=f["ntrans"]
  nplanet=f["nplanet"]
  pair_ttvs=decompose_ttvs(nplanet,ntrans[1:nplanet],pbest_global) .* (24 * 60)
  # println(bestfit)
  n1,n2=ntrans[1],ntrans[2]
  mu1,P1,t01,ecos1,esin1=pbest_global[1:5]
  mu2,P2,t02,ecos2,esin2=pbest_global[6:10]
  # mu3,P3,t03,ecos3,esin3=pbest_global[11:15]
  time1=collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)
  time2=collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  tt1,tt2=tt[1:n1],tt[n1+1:n1+n2]
  tt01,tt02=tt0[1:n1],tt0[n1+1:n1+n2]
  ttmodel1,ttmodel2 = ttmodel[1:n1],ttmodel[n1+1:n1+n2]
  ttsim1,ttsim2=(time1.-t01)./365.25,(time2.-t02)./365.25 #in years
  ttvmodel1,ttvmodel2=(ttmodel1.-time1).*(24*60),(ttmodel2.-time2).*(24*60)
  ttv1,ttv2=(tt1.-time1).* (24 * 60),(tt2.-time2).* (24 * 60) #in minutes
  sigtt1,sigtt2=sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) 
  total1=0;total2=0

  scatter1=(ttvmodel1.-ttv1)
  scatter2=(ttvmodel2.-ttv2)
  println(std(scatter1)," ",std(scatter2))
  # println((tt1.-tt01).* 24*60)

  # xbin,xhist,xbin_square,hist_square=histogram(ttvmodel1.-ttv1,50)
    xbin,xhist,xbin_square,hist_square=histogram([scatter1;scatter2],50)
  plot(xbin_square,hist_square./maximum(hist_square),label=string(fit_type_nplanet))
  # plot(xbin_square,rand(Normal(0, 0.5), n1))
  legend()
end