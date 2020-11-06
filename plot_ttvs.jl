using PyPlot
rc("font",family="serif")
include("decompose_ttvs.jl")

function plot_ttvs(include_moon::Bool=false)
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
    
    figsize=(10,8)
    subplot(211)
    PyPlot.title("Venus TTVs and their sources")
    tick_params(direction="in")
    plot(ttsim1,ttv1,linewidth=1.5,color="grey",label="Total")
    plot(ttsim1,pair_ttvs[1,3,1:n1],linewidth=1.25,color="firebrick",label="Jupiter")
    plot(ttsim1,pair_ttvs[1,2,1:n1],linewidth=1.25,label="Earth")
    errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black")
    ylabel("Venus TTVs [minutes]")
    subplot(212)
    tick_params(direction="in")
    plot(ttsim2,ttv2,linewidth=1.5,color="grey",label="Total")
    plot(ttsim2,pair_ttvs[2,3,1:n2],linewidth=1.25,color="firebrick",label="Jupiter")
    plot(ttsim2,pair_ttvs[2,1,1:n2],linewidth=1.25,color="orange",label="Venus")
    if include_moon
        moon = moon_ttvs(ntrans,pbest_global) .* (24 * 60)
        plot(ttsim2,moon,linewidth=1.25,linestyle="--",color="purple")
        ylabel("Earth TTVs [minutes]")
    else
        ylabel("EMB TTVs [minutes]")
    end
    errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black")
    xlabel("Years Observed [N]")
    tight_layout()

#     legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.0)
#     savefig("name.eps")
end

