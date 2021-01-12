using CALCEPH
using PyPlot
rc("font",family="serif")
include("sim_times.jl")

jd1 = 2.4332825e6
jd2 = 2.4515445e6
jdsize = 1000
sigma = 30.0
t0 = range(jd1,stop=jd2-1,length = jdsize)
nyears = (jd2 - jd1)/365.25 

# pva_sun,pva_venus,pva_earth are the body positions,velocities,and ang. momentum computed for the jd1--jd2 range
# sim_times(jd1::Float64,jd2::Float64,jdsize::Int64,addnoise::Bool=false,sigma::Float64=0.0,EMB::Bool=true,seed::Int=42)
tt1,tt2,n_obs,pva_sun,pva_venus,pva_earth = sim_times(jd1,jd2,jdsize,true,sigma,true)

eph = Ephem("INPUTS/planets.dat") ; prefetch(eph)
options = useNaifId + unitDay + unitAU

function plot_orbits(include_moon::Bool=false)
# pva0,pva1,and pva2 are the body positions,velocities,and ang. momentum computed at the transit times found
pva0 = zeros(9,length(tt2))
pva1 = zeros(9,length(tt2))
pva2 = zeros(9,length(tt2))
for i=1:length(tt2)
pva0[1:9,i] = compute(eph,tt1[i],0.5,10,10,options,2)
pva1[1:9,i] = compute(eph,tt1[i],0.5,2,10,options,2)
if include_moon
  pva_2[1:9,i] = compute(eph,t0[i],0.5,3,10,options,2) 
else
pva2[1:9,i] = compute(eph,tt2[i],0.5,399,10,options,2)
end
end
        
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
