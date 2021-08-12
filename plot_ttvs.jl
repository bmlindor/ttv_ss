using PyPlot,Statistics,JLD2
rc("font",family="sans-serif")
include("decompose_ttvs.jl")
# Create plot of observed TTVs vs model fit
function plot_ttvs(jldfit,include_moon::Bool=false)
  tt,tt0,sigtt,ttmodel = jldfit["tt"],jldfit["tt0"],jldfit["sigtt"],jldfit["ttmodel"]
  pbest_global = jldfit["pbest_global"]
  nplanet,ntrans = jldfit["nplanet"],jldfit["ntrans"]
  pair_ttvs = decompose_ttvs(nplanet,ntrans,pbest_global) .* (24 * 60)
  n1,n2 = ntrans[1],ntrans[2]
  mu1,P1,t01,ecos1,esin1 = pbest_global[1:5]
  mu2,P2,t02,ecos2,esin2 = pbest_global[6:10]
  mu3,P3,t03,ecos3,esin3 = pbest_global[11:15]
  time1 = collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)
  time2 = collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  tt1,tt2 = tt[1:n1],tt[n1+1:n1+n2]
  ttsim1,ttsim2 = (ttmodel[1:n1].-t01)./365.25,(ttmodel[n1+1:n1+n2].-t02)./365.25 #in years
  ttv1,ttv2 = (tt1.-time1).* (24 * 60),(tt2.-time2).* (24 * 60) #in minutes
  sigtt1,sigtt2 = sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes
  println("Simulated with σ= ",mean(sigtt1)*60," second noise")

  if include_moon
  	fig=figure(figsize=(7,7))
    subplot(311)
  else
  	fig=figure(figsize=(7,5))
    subplot(211)
  end
  ax1=gca()
  plot(ttsim1,pair_ttvs[1,3,1:n1]+pair_ttvs[1,2,1:n1],linewidth=2,color="grey")#,label="Total variations")
  plot(ttsim1,pair_ttvs[1,3,1:n1],color="firebrick",label="Jupiter contribution",alpha=0.9)
  plot(ttsim1,pair_ttvs[1,2,1:n1],label="Earth contribution",alpha=0.9,linestyle="--")
  # plot(ttsim1,ttv1,linewidth=1.25,color="grey",label="Obs variations")
  errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")
  ylabel("Venus TTVs [min]")
  xlabel("Time Observed [yrs]")
  ax1.minorticks_on()
  ax1.tick_params(which="both",direction="in",top="true",right="true")
  # ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
  # ax1.tick_params(which="minor",direction="in",top="true",right="true",length=2)
  ax1.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
      ncol=3,mode="expand",borderaxespad=0.0)

  if include_moon
    subplot(312,sharex=ax1)
  else 
    subplot(212,sharex=ax1)
  end
  ax2=gca()
  plot(ttsim2,pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2],linewidth=2,color="grey")#,label="Total variations")
  plot(ttsim2,pair_ttvs[2,3,1:n2],color="firebrick",label="Jupiter contribution",alpha=0.9)
  plot(ttsim2,pair_ttvs[2,1,1:n2],color="orange",label="Venus contribution",alpha=0.9,linestyle="--")
  ylabel("EMB TTVs [min]")
  errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")
  xlabel("Time Observed [yrs]")
  ax2.minorticks_on()
  ax2.tick_params(which="both",direction="in",top="true",right="true")
  ax2.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
      ncol=3,mode="expand",borderaxespad=0.0)

  if include_moon
    jldfit2=jldopen("FITS/moon_fit30.0s30.0yrs.jld2","r")
    tt,tt0,sigtt,ttmodel = jldfit2["tt"],jldfit2["tt0"],jldfit2["sigtt"],jldfit2["ttmodel"]
    pbest_global = jldfit2["pbest_global"]
    nplanet,ntrans = jldfit2["nplanet"],jldfit2["ntrans"]
    pair_ttvs = decompose_ttvs(nplanet,ntrans,pbest_global) .* (24 * 60)
    n1,n2 = ntrans[1],ntrans[2]
    mu1,P1,t01,ecos1,esin1 = pbest_global[1:5]
    mu2,P2,t02,ecos2,esin2 = pbest_global[6:10]
    mu3,P3,t03,ecos3,esin3 = pbest_global[11:15]
    time1 = collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)
    time2 = collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
    tt1,tt2 = tt[1:n1],tt[n1+1:n1+n2]
    ttsim1,ttsim2 = (ttmodel[1:n1].-t01)./365.25,(ttmodel[n1+1:n1+n2].-t02)./365.25 #in years
    ttv1,ttv2 = (tt1.-time1).* (24 * 60),(tt2.-time2).* (24 * 60) #in minutes
    sigtt1,sigtt2 = sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes
  	moon = moon_ttvs(ntrans,pbest_global) .* (24 * 60)
    subplot(313,sharex=ax1)
    ax3=gca()
    plot(ttsim2,pair_ttvs[2,3,1:n2]+pair_ttvs[2,1,1:n2]+moon,linewidth=2,color="grey")#,label="Total variations")
    plot(ttsim2,pair_ttvs[2,3,1:n2],color="firebrick",label="Jupiter contribution",alpha=0.9)
    plot(ttsim2,pair_ttvs[2,1,1:n2],color="orange",label="Venus contribution",alpha=0.9,linestyle="--")
    plot(ttsim2,moon,linestyle="-.",color="purple",label="Moon contribution",alpha=0.9)
    ylabel("Earth TTVs [min]")
    errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",mec="black",mfc="white")
    xlabel("Time Observed [yrs]")
    ax3.minorticks_on()
    ax3.tick_params(which="both",direction="in",top="true",right="true")
    # ax3.legend(loc="lower right")
  	ax3.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc="lower left",
        ncol=3,mode="expand",borderaxespad=0.0)
  end
  tight_layout()
  savefig("IMAGES/ttvs.eps")
  show()
end