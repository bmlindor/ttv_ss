#Julia v1.1
using PyPlot, CALCEPH, DelimitedFiles
using Statistics, DataFitting, Random
using Unitful, UnitfulAstro, LinearAlgebra
if !@isdefined(CGS)
  include("CGS.jl")
  using Main.CGS
end
include("regress.jl")

# Initial JD times for days in 100 years
nyear = 100
np0 = 365*nyear 
t0 = 2451544.5 - 50*365.25 .+ range(0.5,stop = np0 - 0.5,length = np0)

# Make a circle to represent the Sun:
theta_sun = range(0,stop = 2*pi,length = 100)
xsun = CGS.RSUN/CGS.AU * cos.(theta_sun)
ysun = CGS.RSUN/CGS.AU * sin.(theta_sun)

# For planets without moons, Mercury and Venus, 
# the barycenter location coincides with the body center of mass. 

# CALCELPH retrieves the position, velocity and acceleration of Earth (geocenter) relative
# to the Earth-Moon system barycenter in kilometers, kilometers per second and
# kilometers per second square at JD= 2451624.5 TDB timescale
# for best accuracy the first time argument should be the integer part and the
# delta the fractional part 
# adjust as step through day? 

eph = Ephem("planets.dat") ; prefetch(eph)
options = useNaifId + unitDay + unitAU

pva_sun = zeros(9, np0)
pva_venus = zeros(9, np0)
pva_earth = zeros(9, np0)

for i=1:np0
   pva_sun[1:9,i] = compute(eph,t0[i],0.5,10,10,options,2)
   pva_venus[1:9,i] = compute(eph,t0[i],0.5,2,10,options,2) # useNaifId = 2 for Venus 
   pva_earth[1:9,i] = compute(eph,t0[i],0.5,3,10,options,2) # useNaifId = 3 for Earth
end

L_earth = cross(pva_earth[1:3], pva_earth[4:6])
L_venus = cross(pva_venus[1:3], pva_venus[4:6])

n_obs = cross(L_earth,L_venus)
n_obs /= norm(n_obs) #from one direction when both transit
x_obs, y_obs, z_obs = n_obs[1], n_obs[2], n_obs[3]

# Finds the transit by calculating the position, velocity, and acceleration for a body,
    # Finds local minimum of f(t) by solving for time when f(t) = x_bar dot v_bar = 0
    # Returns transit time in JD
function find_transit(body_id, eph, jd1, jd2, n_obs, N)
    # N = jd2 - jd1
    JD_0 = 0.0
    ff = zeros(N)
    xdotn = 0.0
    pos = zeros(3,N) # position of body relative to Sun
    # Compute functions of position and velocity wrt time:
    function calc_ffs(t)
        pva = compute(eph,JD_0,t,body_id,10,options,2)
        x = pva[1:3]; v = pva[4:6]; a = pva[7:9];
        f = dot(x, v) - (dot(x, n_obs))*(dot(v, n_obs))
        Df = dot(v, v) + dot(x, a) - (dot(v, n_obs))^2 - (dot(x, n_obs)*dot(a, n_obs))
        #Df *= ( 3600 * 24 / 1) # Converts to units of days
        return f, Df, dot(x, n_obs), x
    end
    # Computing minimum sky separation of planet wrt star for all JDs
    dt =  (jd2 - jd1)/(N-1)
    JD = zeros(N)
    i_min = 1
    ff_min = Inf
    for i=1:N
        JD[i] = jd1 + dt*(i-1)
        JD_0 = JD[i]
        ff[i], Df, xdotn, pos[:,i] = calc_ffs(0.0)
    # Estimate of transit time:
            # Df > 0 for transit occuring; 
            # xdotn > 0 for planet in front of star as seen by observer; 
            # local minimum value over entire range (i.e. when close to zero)
        if (Df > 0) && (xdotn > 0) && (abs(ff[i]) < ff_min)
            i_min = i 
            ff_min = abs(ff[i])
        end

    end
     #println("Estimated Transit Time: ", JD[i_min])
    # Refine initial guess using linear approx:
    JD_0 = JD[i_min]
    JD_n = 0.0

    JD_n1 = JD_n + 1
    JD_n2 = JD_n + 2
    iter = 0
    # ITMAX = 20
    ITMAX = 6 
    # we've found that we don't need large ITMAX to find solution; does that change for different bodies?
    for iter=0:ITMAX
        JD_n2 = JD_n1
        JD_n1 = JD_n
        while JD_n > 1
            JD_n -= 1.0
            JD_0 += 1.0
        end
        while JD_n < 0
            JD_n += 1.0
            JD_0 -= 1.0

        end
        f_n, Df_n, xdotn, x = calc_ffs(JD_n)
        JD_n -= f_n/Df_n 

        # Break out if we have reached maximum iterations, or if
              # current transit time estimate equals one of the prior two steps:
        if (JD_n == JD_n1) || (JD_n == JD_n2)
            break
        end
    end          

    JD_tt = JD_0 + JD_n
    #println("Refined Transit Time: ", JD_tt)
    return JD, ff, i_min, pos, JD_tt
end

# Find the transit times for given body_id, planetary period, and number of refinement steps N
function find_times(body_id, eph, t0, period, period_err, n_obs, N)
    times = Float64[]
    t_final = t0[end]
    i=1
    # initializes & finds first transit time
    JD,ff,i_min,pos,JD_tt = find_transit(body_id,eph,t0[i],t0[i]+period,n_obs,1000)
    push!(times, JD_tt)
    
    # Find subsequent transit times by shifting time frame by 1 planetary period
    while JD_tt < t_final
        t_start = JD_tt+period-period_err
        t_end = JD_tt+period+period_err
        JD,ff,i_min,pos,JD_tt = find_transit(body_id,eph,t_start,t_end,n_obs,N)
        push!(times, JD_tt)
    end
    return times
end

P_earth = 365
P_venus = 225
P_err = 2

tt_earth = find_times(3, eph, t0, P_earth, P_err, n_obs, 10)
tt_venus = find_times(2, eph, t0, P_venus, P_err, n_obs, 10)

ntearth = length(tt_earth)
ntvenus = length(tt_venus)

# Find ttvs via linear regression of transit time data
# accounts for missing transits (noncontinuous) 
# by rounding [difference in consecutive transit times/Period]
function find_ttvs(tt, period; sigma_x = ones(length(tt)))
    nt = length(tt)
    x = zeros(2,nt)
    x[1,1:nt] .= 1.0
    x[2,1] = 0.0 # for fitting time of first transit
    for i=2:nt
        x[2,i] = x[2,i-1] + round((tt[i]-tt[i-1])/period) 
    end
        
    coeff, cov = regress(x, tt, sigma_x)
    # coeff[1] = best linear fit approx of first tt, coeff[2] = average period
    ttv = tt .- coeff[1].*vec(x[1,1:nt]) .- coeff[2].*vec(x[2,1:nt])
    return coeff, ttv
end

coeff_earth, ttv_earth = find_ttvs(tt_earth, P_earth)
coeff_venus, ttv_venus = find_ttvs(tt_venus, P_venus)

writedlm("ttv_earth.txt", zip(tt_earth,ttv_earth))
writedlm("ttv_venus.txt", zip(tt_venus,ttv_venus))
# need to include ttv wrappers
# function chisquare2(param)
#     chisq = 0.0
#     tt_model = ttv_wrapper2(tt0,param)
#     for j=1:length(tt0)
#       chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
#     end
#     return chisq
# end

# function chisquare3(param)
#     chisq = 0.0
#     tt_model = ttv_wrapper3(tt0,param)
#     for j=1:length(tt0)
#       chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
#     end
#     return chisq
# end

function plot_ecliptic()
    title("Orbits Along Ecliptic")
    plot(xsun,ysun,label="Sun")
    plot(vec(pva_venus[2,1:np0]),vec(pva_venus[3,1:np0]),label="Venus")
    plot(vec(pva_earth[2,1:np0]),vec(pva_earth[3,1:np0]),label="Earth")
    xlabel("[AU]")
    ylabel("[AU]")
    legend()
end

# function plot_ttvs()
# end


# Okay, so now add noise to the TTVs of both bodies:

# E_ttv = readdlm("agol_ttv_earth.txt") #Earth
# V_ttv = readdlm("agol_ttv_venus.txt") #Venus
# #Transit Times
# E_tt = E_ttv[:,1] .* 24*60; V_tt = V_ttv[:,1].* 24*60; #days --> minutes
# #Transit Timing Variations
# E_tv = E_ttv[:,2] .* 24 * 60 ; V_tv = V_ttv[:,2] .*24 * 60 ;#days --> minutes
# E_tt

# dom1 = Domain(tt1) ; dom2 = Domain(V_tt)
# rng = MersenneTwister(0)
# E_nvar = Statistics.std(E_tv) ; V_nvar = Statistics.std(V_tv)
# E_noise = randn(rng, length(E_tv)) * E_nvar/2 ; V_noise = randn(rng, length(V_tv)) * V_nvar/2
# E_tv_noised = E_noise + E_tv ; V_tv_noised = V_noise + V_tv
# data1 = Measures(E_tv_noised, 1.0) ; data2 = Measures(V_tv_noised, 1.0)

# subplot(211) 
# suptitle("Transit Timing Variations") # Supe title, title for all subplots combine
# plot(V_tt, V_tv, color = "k")
# plot(V_tt, V_tv_noised, label = "+noise")
# ylabel("Venus [min]"); legend()
# subplot(212)
# plot(E_tt, E_tv, color = "k")
# plot(E_tt, E_tv_noised, label = "+noise")
# ylabel("Earth [min]"); legend()
# xlabel("transit times")
# #ylabel("timing variations [minutes]") #legend();
# plot(E_tt, E_tv, label = "Earth", color = "k")
# errorbar(E_tt,E_tv,yerr = E_nvar, label ="Std", fmt = ".", alpha = 0.75)
# legend();xlabel("transit times");ylabel("timing variations [minutes]")
# data1 = readdlm("ttv_venus.txt")
# tt1 = vec(data1[:,1]); ttv1 = vec(data1[:,2])

# data2 = readdlm("ttv_earth.txt")
# tt2 = vec(data2[:,1]); ttv2 = vec(data2[:,2])

# nt1 = length(tt1)
# t01 = coeff_venus[1]; per1 = coeff_venus[2]

# nt2 = length(tt2)
# t02 = coeff_earth[1]; per2 = coeff_earth[2];
# # best fit linear transit times w/o ttvs
# t1  = collect(t01 .+ per1 .* range(0,stop = nt1-1,length = nt1)) 
# t2  = collect(t02 .+ per2 .* range(0,stop = nt2-1,length = nt2))
#     # Best-fit linear transit times:
# tt0 = [t1;t2] # appends t2 times to t1 times
# weight = ones(nt1+nt2) #assigns each data point stat weight d.t. noise [will add noise later]
# # Actual transit times:
# tt = [tt1;tt2]

# Okay, now let's do a 2-planet fit:
# initial params defined for both planets
# mass ratio, period, initial transit time, e*cos(omega), e*sin(omega)
# param = [3e-6,per1,t01,0.01,0.01,3e-6,per2,t02,0.01,0.01] 
# println("Initial parameters: ",param)

# println("Need to define ntrans and nplanet!!")
# nplanet = 2
# ntrans = [nt1,nt2]

