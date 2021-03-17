using PyPlot
rc("font",family="sans-serif")
include("decompose_ttvs.jl")


function moon_contribution()
  fig, ax1 = subplots(figsize=(8,4))
  plot(ttsim2,moon,linewidth=1.25,linestyle="--",color="purple",label="Moon's contribution")
  ylabel("Earth TTVs [minutes]")
  xlabel("Years Observed [N]")
  ax1.minorticks_on()
  ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
  ax1.tick_params(which="minor",direction="in",top="true",right="true",length=2)
  tight_layout()
  show()
end

function plot_ttvs(jldfit,include_moon::Bool=false)
  tt,tt0,sigtt,ttmodel = jldfit["tt"],jldfit["tt0"],jldfit["sigtt"],jldfit["ttmodel"]
  pbest_global = jldfit["pbest_global"]
  nplanet,ntrans = jldfit["nplanet"],jldfit["ntrans"]
  pair_ttvs = decompose_ttvs(nplanet,ntrans,pbest_global) .* (24 * 60)
  n1,n2,n3 = ntrans
  mu1,P1,t01,ecos1,esin1 = pbest_global[1:5]
  mu2,P2,t02,ecos2,esin2 = pbest_global[6:10]
  mu3,P3,t03,ecos3,esin3 = pbest_global[11:15]
  time1 = collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)
  time2 = collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  tt1 = tt[1:n1]
  tt2 = tt[n1+1:n1+n2]
  ttsim1 = (ttmodel[1:n1].-t01)./365.25 #in years
  ttsim2 = (ttmodel[n1+1:n1+n2].-t02)./365.25 
  ttv1 = (tt1.-time1).* (24 * 60) #in minutes
  ttv2 = (tt2.-time2).* (24 * 60) 
  sigtt1 = sigtt[1:n1].* (24 * 60)
  sigtt2 = sigtt[n1+1:n1+n2].* (24 * 60)
  moon = moon_ttvs(ntrans,pbest_global) .* (24 * 60)
  # fig, ax1 = subplots(figsize=(8,3))
  fig=figure(figsize=(8,6))
  subplot(211)
  ax1=gca()
  plot(ttsim1,pair_ttvs[1,3,1:n1],linewidth=1.25,color="firebrick",label="Jupiter's contribution")
  plot(ttsim1,pair_ttvs[1,2,1:n1],linewidth=1.25,label="Earth's contribution")
  plot(ttsim1,ttv1,linewidth=1.5,color="grey",label="Total TTVs")
  errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black")
  ylabel("Venus TTVs [minutes]")
  xlabel("Years Observed [N]")
  ax1.minorticks_on()
  ax1.tick_params(which="both",direction="in",top="true",right="true")
  # ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
  # ax1.tick_params(which="minor",direction="in",top="true",right="true",length=2)
  ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left",
         ncol=2, mode="expand", borderaxespad=0.0)
  # tight_layout()
  # show()
  # savefig("IMAGES/venusttvs.png")

  # fig, ax2 = subplots(figsize=(8,3))
  subplot(212,sharex=ax1)
  ax2=gca()
  plot(ttsim2,pair_ttvs[2,3,1:n2],linewidth=1.25,color="firebrick",label="Jupiter's contribution")
  plot(ttsim2,pair_ttvs[2,1,1:n2],linewidth=1.25,color="orange",label="Venus's contribution")
  if include_moon
    plot(ttsim2,moon,linewidth=1.25,linestyle="--",color="purple",label="Moon's contribution")
    ylabel("Earth TTVs [minutes]")
  else
    ylabel("EMB TTVs [minutes]")
  end
  plot(ttsim2,ttv2,linewidth=1.5,color="grey",label="Total TTVs")
  errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black")
  xlabel("Years Observed [N]")
  ax2.minorticks_on()
  ax2.tick_params(which="major",direction="in",top="true",right="true",length=6)
  ax2.tick_params(which="minor",direction="in",top="true",right="true",length=2)
  ax2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left",
     ncol=2, mode="expand", borderaxespad=0.0)
  tight_layout()
  show()
#     legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.0)
end



