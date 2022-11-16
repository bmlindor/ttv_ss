using Test,CALCEPH,PyPlot
include("sim_times.jl")
eph = Ephem("/Users/bethleelindor/work/washington/premap2022/jup344.bsp") ; prefetch(eph)
options = useNaifId+unitKM+unitDay # useNaifId + unitDay + unitAU
AU = 149597870.700 #km
jd1=2.4332825e6
jdsize=10
function test_sim_obs_and_find_times()
  sigma=30
  nyear=10
  obs="fromEMB"
  body,tt0,tt,sigtt=sim_obs_and_find_times(jd1,sigma,nyear,obs)
end
# body,tt0,tt,sigtt=test_sim_obs(jd1,sigma,nyear,obs)

# JD_venus,ff_venus,i_min_venus,pos_venus,tt_venus = find_transit(2,eph,jd1,jd1+365,n_obs,365)
# JD_earth,ff_earth,i_min_earth,pos_earth,tt_earth = find_transit(3,eph,jd1,jd1+365,n_obs,365)

# P_venus = 225.0
# P_earth = 365.0
# P_err = 1.0
# tt1 = transit_times(2,eph,t0,P_venus,P_err,n_obs,10)
# tt2 = transit_times(3,eph,t0,P_earth,P_err,n_obs,10)
# subplot(211)
# plot(JD_venus,ff_venus)
# plot(JD_earth,ff_earth)
# scatter(tt_venus,tt_venus)
# scatter(tt_earth,tt_earth)
# subplot(212)
# plot(,tt_earth)
# plot(,tt_venus)
# show()
# @testset "TTV_SS" begin
#     @test test_sim_obs_and_find_times()
# end