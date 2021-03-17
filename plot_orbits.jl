using CALCEPH,PyPlot
rc("font",family="sans-serif")
include("sim_times.jl")

function plot_orbits(include_moon::Bool=false)

jd1 = 2.4332825e6
nyear = 40.0
jd2 = nyear*365.25 + jd1
jdsize = 1000
sigma = 30.0
t0 = range(jd1,stop=jd2-1,length = jdsize)

# pva_sun,pva_venus,pva_earth are the body positions,velocities,and ang. momentum computed for the jd1--jd2 range
#sim_times(jd1::Float64,nyear::Float64,addnoise::Bool=false,sigma::Float64=0.0,EMB::Bool=true,seed::Int=42)
tt1,tt2,n_obs,pva_sun,pva_venus,pva_earth = sim_times(jd1,nyear,true,sigma,true)

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
	if include_moon
	  pva_2[1:9,i] = compute(eph,t0[i],0.5,3,10,options,2)./AU
	else
		pva2[1:9,i] = compute(eph,tt2[i],0.5,399,10,options,2)./AU
	end
end

function 2D_orbits()
fig = plt.figure(figsize=(6, 6))
ax1 = fig.add_subplot(211)
ax1.plot(vec(pva_sun[2,:]),vec(pva_sun[3,:]),label="Sun",color=:yellow,marker="o")
# ax1.plot(vec(pva_sun[2,:]),vec(pva_sun[3,:]),color=:yellow,marker="o",mec="black")
ax1.plot(vec(pva_venus[2,:]),vec(pva_venus[3,:]),label="Venus",color=:orange)
ax1.plot(vec(pva_earth[2,:]),vec(pva_earth[3,:]),label="Earth")
ax1.plot([0,n_obs[2]*1.1],[0,n_obs[3]*1.1],color=:black)
ax1.tick_params(which="major",direction="in",length=6,
    left="false",right="false",top="false",bottom="false",
    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
xlabel("[AU]")
ylabel("[AU]")
# legend(loc="lower left")

subplot(212,sharex=ax1)
ax2=gca()
ax2.scatter(vec(pva0[2,:]),vec(pva0[3,:]),color=:yellow,marker="o",edgecolors="black")
# ax2.scatter(pva_venus[1,imin1],pva_venus[2,imin1],label="Venus Transit",color=:orange)
ax2.plot(vec(pva_venus[1,:]),vec(pva_venus[2,:]),color=:gray)
ax2.scatter(vec(pva1[1,:]),vec(pva1[2,:]),label="Venus",color=:orange,marker=".")
# ax2.scatter(pva_earth[1,imin2],pva_earth[2,imin2],label="Earth Transit")
ax2.plot(vec(pva_earth[1,:]),vec(pva_earth[2,:]),color=:gray)
ax2.scatter(vec(pva2[1,:]),vec(pva2[2,:]),label="Earth",marker=".")
ax2.plot([0,n_obs[1]*1.1],[0,n_obs[2]*1.1],"k--")
ax2.tick_params(which="major",direction="in",
    left="true",right="false",top="false",bottom="true",
    labelbottom="true",labeltop="false",labelleft="true",labelright="false")
ax2.legend()
xlabel("[AU]")
ylabel("[AU]")
# savefig("sim_times.eps")
end

function 3D_orbits()
fig=figure(figsize=(8,6))
# PyPlot.scatter3D(1,0,0,color=:black,marker="s")
# PyPlot.text3D(1.1,0,0,"x")
# PyPlot.scatter3D(1,0,0,color=:black,marker=">")
# PyPlot.scatter3D(0,-1.2,0,color=:black,marker="<")
# PyPlot.scatter3D(0,0,1,color=:black,marker="^")
# PyPlot.plot3D([0,1*1.5],[0,0*1.5],[0,0*1.5],linewidth=3,color=:gray)
# PyPlot.plot3D([0,0*1.5],[0,-1*1.5],[0,0*1.5],linewidth=3,color=:gray)
# PyPlot.plot3D([0,0*1.5],[0,0*1.5],[0,1*1.5],linewidth=3,color=:gray)
# PyPlot.scatter3D(xsun,ysun,0,marker="o",color=:yellow,ms=20)
PyPlot.plot3D(vec(pva_venus[1,:]), vec(pva_venus[2,:]), vec(pva_venus[3,:]),color=:orange,label="Venus Orbit",alpha=0.25)
PyPlot.plot3D(vec(pva_earth[1,:]), vec(pva_earth[2,:]), vec(pva_earth[3,:]),color=:skyblue,label="Earth Orbit",alpha=0.5)
PyPlot.plot3D([0,n_obs[1]*1.2],[0,n_obs[2]*1.2],[0,n_obs[3]*1.2],linestyle="--",color=:black,alpha=0.5)
PyPlot.plot3D(vec(pva0[1,:]), vec(pva0[2,:]), vec(pva0[3,:]),color=:yellow,marker="o",ms=20,mec=:gold)
PyPlot.plot3D(vec(pva1[1,:]),vec(pva1[2,:]),vec(pva1[3,:]),color=:orange,marker="o")
PyPlot.plot3D(vec(pva2[1,:]),vec(pva2[2,:]),vec(pva2[3,:]),marker="o")
# ax.tick_params(which="minor",direction="in",length=2,
#     left="false",right="false",top="true",bottom="true",
#     labelbottom="false",labeltop="false",labelleft="false",labelright="false")
# xlim(-1,1)
# ylim(-1,1)
xlabel("x [AU]")
ylabel("y [AU]")
zlabel("z [AU]")    
end
end
