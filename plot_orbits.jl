using CALCEPH
include("sim_times.jl")

jd1 = 2.4332825e6
jd2 = 2.4515445e6
jdsize = 1000
sigma = 30.0
t0 = range(jd1,stop=jd2-1,length = jdsize)
nyears = (jd2 - jd1)/365.25 

pva1, pva2, pva0, imin1, imin2, n_obs = sim_times(jd1,jd2,jdsize,true,sigma,true)

function plot_orbits()
subplot(211)
title("Orbits Along Ecliptic")
plot(xsun,ysun,label="Sun")
plot(vec(pva_venus[2,1:jdsize]),vec(pva_venus[3,1:jdsize]),label="Venus")
plot(vec(pva_earth[2,1:jdsize]),vec(pva_earth[3,1:jdsize]),label="Earth")
xlabel("[AU]")
ylabel("[AU]")
legend()
subplot(212)
title("Top-Down Orbits")
plot(xsun,ysun,label="Sun")
plot(vec(pva_venus[1,1:jdsize]),vec(pva_venus[2,1:jdsize]),label="Venus")
plot(vec(pva_earth[1,1:jdsize]),vec(pva_earth[2,1:jdsize]),label="Earth")
xlabel("[AU]")
ylabel("[AU]")
legend(loc="lower left")
clf()

test = 365 
i=1
JD_venus,ff_venus,i_min_venus,pos_venus,tt_venus = find_transit(2,eph,t0[i],t0[i]+test,n_obs,test)
JD_earth,ff_earth,i_min_earth,pos_earth,tt_earth = find_transit(3,eph,t0[i],t0[i]+test,n_obs,test)
# title("Top-Down Orbits w/ Observer")
figsize=(8,8)
plot(pos_venus[1,i_min_venus],pos_venus[2,i_min_venus],"o",label="Venus Transit",color=:orange)
plot(pos_earth[1,i_min_earth],pos_earth[2,i_min_earth],"o",label="Earth Transit")
plot(xsun,ysun,"o",label="Sun",color=:yellow)
plot(pos_venus[1,:],pos_venus[2,:],color=:grey)
plot(pos_earth[1,:],pos_earth[2,:],color=:grey)
plot([0,x_obs*1.1],[0,y_obs*1.1],"k--")
legend(loc="upper left")
xlabel("[AU]")
ylabel("[AU]")
# savefig("sim_times.eps")
end
