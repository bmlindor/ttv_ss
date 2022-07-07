using CALCEPH,PyPlot
rc("font",family="sans-serif")
include("sim_times.jl")

function plot_orbits(jd1::Float64,sigma::Real,nyear::Real,dimension::Int64,EMB::Bool=false)

jd1 = 2.4332825e6
# nyear = 40.0
jd2 = nyear*365.25 + jd1
jdsize = 1000
# sigma = 30.0
t0 = range(jd1,stop=jd2-1,length = jdsize)

# pva_sun,pva_venus,pva_earth are the body positions,velocities,and ang. momentum computed for the jd1--jd2 range
#sim_times(jd1::Float64,nyear::Float64,addnoise::Bool=false,sigma::Float64=0.0,EMB::Bool=true,seed::Int=42)
tt1,tt2,n_obs,pva_sun,pva_venus,pva_earth = sim_times(jd1,sigma,nyear,EMB)

eph = Ephem("INPUTS/DE440.bsp") ; prefetch(eph)
options = useNaifId+unitKM+unitDay # useNaifId + unitDay + unitAU
AU = 149597870.700 #km

theta_sun = range(0,stop = 2pi,length = 100)
xsun = CGS.RSUN/CGS.AU * cos.(theta_sun)
ysun = CGS.RSUN/CGS.AU * sin.(theta_sun)
# pva0,pva1,and pva2 are the body positions,velocities,and ang. momentum computed at the transit times found
pva0 = zeros(9,length(tt2))
pva1 = zeros(9,length(tt2))
pva2 = zeros(9,length(tt2))
for i=1:length(tt2)
pva0[1:9,i] = compute(eph,tt1[i],0.5,10,10,options,2)./AU
pva1[1:9,i] = compute(eph,tt1[i],0.5,2,10,options,2)./AU
	if EMB
	  pva_2[1:9,i] = compute(eph,t0[i],0.5,3,10,options,2)./AU
	else
		pva2[1:9,i] = compute(eph,tt2[i],0.5,399,10,options,2)./AU
	end
end

if dimension==2
fig = plt.figure(figsize=(6, 6))
ax1 = fig.add_subplot(111)
# ax1.plot(vec(pva_sun[2,:]),vec(pva_sun[3,:]),color=:yellow,marker="o",mec="black")
ax1.plot(vec(pva_venus[1,:]),vec(pva_venus[2,:]),alpha=0.25)
ax1.plot(vec(pva_earth[1,:]),vec(pva_earth[2,:]),alpha=0.25)
ax1.plot(vec(pva1[1,:]),vec(pva1[2,:]),label="Venus Transit",marker="o")
ax1.plot(vec(pva2[1,:]),vec(pva2[2,:]),color=:orange,label="Earth Transit",marker="o")
ax1.plot([0,n_obs[1]*1.1],[0,n_obs[2]*1.1],color=:grey,linestyle="--",alpha=0.5)
ax1.plot(xsun,ysun,color=:yellow,marker="o",ms=5,mec=:gold)
# ax1.plot(vec(pva_sun[2,:]),vec(pva_sun[3,:]),label="Sun",color=:yellow,marker="o",ms=10,mec="gold")
ax1.tick_params(which="major",direction="in",length=6,
    left="false",right="false",top="false",bottom="false",
    labelbottom="true",labeltop="false",labelleft="true",labelright="false")
# legend()
xlabel("x-position [AU]")
ylabel("y-position [AU]")
# legend(loc="lower left")
end

if dimension==3
fig=figure(figsize=(8,6))
# PyPlot.scatter3D(xsun,ysun,0,marker="o",color=:yellow,ms=20)
PyPlot.plot3D(vec(pva_venus[1,:]), vec(pva_venus[2,:]), vec(pva_venus[3,:]),alpha=0.25)
PyPlot.plot3D(vec(pva_earth[1,:]), vec(pva_earth[2,:]), vec(pva_earth[3,:]),alpha=0.25)
PyPlot.plot3D([0,n_obs[1]*1.2],[0,n_obs[2]*1.2],[0,n_obs[3]*1.2],linestyle="--",color=:grey,alpha=0.5)
PyPlot.plot3D(vec(pva0[1,:]), vec(pva0[2,:]), vec(pva0[3,:]),color=:yellow,marker="o",ms=10,mec=:gold)
PyPlot.plot3D(vec(pva1[1,:]),vec(pva1[2,:]),vec(pva1[3,:]),label="Venus",marker="o")
PyPlot.plot3D(vec(pva2[1,:]),vec(pva2[2,:]),vec(pva2[3,:]),color=:orange,label="Earth",marker="o")
PyPlot.tick_params(which="major",
    left="false",right="false",top="false",bottom="false")
legend()
# xlim(-1,1)
# ylim(-1,1)
xlabel("x [AU]")
ylabel("y [AU]")
zlabel("z [AU]")    
end
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