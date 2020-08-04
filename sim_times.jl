# Julia v1.1
using PyPlot, CALCEPH 
using FITSIO, DelimitedFiles
using Statistics, DataFitting, Random
using Unitful, UnitfulAstro, LinearAlgebra
if !@isdefined(CGS)
  include("CGS.jl")
  using Main.CGS
end
include("regress.jl")

function sim_times(jd1::Float64, jd2::Float64, Nsteps::Int64, addnoise::Bool=false, sigma::Float64=0.0, seed::Int=42)
    # To do: output file with arguments in header
    #       output data frame with header?

    # Initial JD times for days in 100 years
    # nyear = 100
    # np0 = 365*nyear 
    @assert jd1 >= 2414105.0
    @assert jd2 <= 2488985.0
    Random.seed!(seed)
    dt = (jd2 - jd1)/Nsteps
    t0 = range(jd1, stop=jd2-1, length = Nsteps)
    # t0 = 2451544.5 - 50*365.25 .+ range(0.5,stop = np0 - 0.5,length = np0)
    # println(t0[1]) #= 2.4332825e6  

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

    # Load ephemerides from data and set units
    eph = Ephem("planets.dat") ; prefetch(eph)
    options = useNaifId + unitDay + unitAU

    # Find observer location required to see transits
    # pva_sun = zeros(9, np0)

    pva_sun = zeros(9, Nsteps)
    pva_venus = zeros(9, Nsteps)
    pva_earth = zeros(9, Nsteps)
    # for i=1:np0
    for i=1:Nsteps
       pva_sun[1:9,i] = compute(eph,t0[i],0.5,10,10,options,2)
       pva_venus[1:9,i] = compute(eph,t0[i],0.5,2,10,options,2) # useNaifId = 2 for Venus 
       pva_earth[1:9,i] = compute(eph,t0[i],0.5,3,10,options,2) # useNaifId = 3 for Earth-Moon bary
    end
    L_venus = cross(pva_venus[1:3], pva_venus[4:6])
    L_earth = cross(pva_earth[1:3], pva_earth[4:6])
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
    P_venus = 225
    P_earth = 365
    P_err = 2
    tt1 = find_times(2, eph, t0, P_venus, P_err, n_obs, 10)
    nt1 = length(tt1)
    tt2 = find_times(3, eph, t0, P_earth, P_err, n_obs, 10)
    nt2 = length(tt2)
    # Actual transit times:

    # Find ttvs via linear regression of transit time data
    # accounts for missing transits (noncontinuous) 
    # by rounding [difference in consecutive transit times/Period]
    # function find_ttvs(tt, period; sigma_x = ones(length(tt)))
    function find_coeffs(tt, period, addnoise, sigma)
        nt = length(tt)
        noise = zeros(nt)
        x = zeros(2,nt)
        x[1,1:nt] .= 1.0
        x[2,1] = 0.0 # for fitting time of first transit
        for i=2:nt
            x[2,i] = x[2,i-1] + round((tt[i]-tt[i-1])/period) 
        end
        # coeff, cov = regress(x, tt, sigma_x)
        # Add noise to transit times
        if addnoise
            sigtt = ones(nt) * sigma / (24 * 3600) # sigma in seconds, sigtt in days
            noise = randn(nt) .* sigtt  
            # println("Noise added with σ of ", string(sigma), " seconds.")
        else
            sigtt = ones(nt)
            # println("No noise added.")
        end
        coeff, covcoeff = regress(x, tt, sigtt)
        # println(tt, sigtt, std(sigtt))
        # coeff[1] is best linear fit approx of first tt, coeff[2] is average period
        ttv = tt .- coeff[1].*vec(x[1,1:nt]) .- coeff[2].*vec(x[2,1:nt])
        # return coeff, ttv
        # might want to return sigtt at some point?
        return coeff, noise, sigtt,ttv
    end
    # coeff_venus, ttv1 = find_ttvs(tt1, P_venus,noise::Bool)
    coeff_venus, noise1, sigtt1, ttv1 = find_coeffs(tt1, P_venus, addnoise, sigma);
    coeff_earth, noise2, sigtt2, ttv2 = find_coeffs(tt2, P_earth, addnoise, sigma);

    t01 = coeff_venus[1]; per1 = coeff_venus[2]
    t02 = coeff_earth[1]; per2 = coeff_earth[2]

    # best fit linear transit times w/o ttvs
    t1  = collect(t01 .+ per1 .* range(0,stop = nt1-1,length = nt1)) 
    t2  = collect(t02 .+ per2 .* range(0,stop = nt2-1,length = nt2))
    # Best-fit linear transit times:
    tt0 = [t1;t2] # appends t2 times to t1 times


    # Plot transit times and TTVs
    function plot_ttvs(sigma)
        # subplot(211)
        # scatter((t1.-t01)./per1,tt1.-t1) #x is tranit number 
        # plot((t1.- t01)./per1,ttv1) 
        # errorbar((t1.-t01)./per1,ttv1, noise1)
        scatter((t1.-t01)./365.25, tt1.-t1) # x is JD in years
        plot((t1.-t01)./365.25, ttv1)
        # subplot(212)
        scatter((t2.-t02)./365.25, tt2.-t2, color="green")
        plot((t2.-t02)./365.25, ttv2)
        errorbar((t2.-t02)./365.25,ttv2, noise2)
        # scatter((t2.-t02)./per2,tt2.-t2,color="green") 
        # plot((t2.-t02)./per2,ttv2)
        # title(sigma)
        xlabel("JD (years)")
        ylabel("TTVs")
        # savefig("OUTPUTs/")
    end
    # plot_ttvs(sigma)

    function write_file(jd1, jd2, sigma)
        body = zeros((nt1+nt2))
        body[1:nt1] .= 1.0
        body[nt1+1:nt1+nt2] .= 2.0
        tt = [tt1+noise1;tt2+noise2]
        if addnoise
            noise = [noise1;noise2]
            sigtt = [sigtt1;sigtt2]
            name = string("INPUTS/tt_data",sigma,"s.txt")
            writedlm(name, zip(body, tt0, tt, sigtt))
        else
            writedlm("INPUTS/tt_data.txt", zip(body, tt))
        end
    end
    write_file(jd1, jd2, sigma)
    return tt1, noise1+tt1, tt2, noise2+tt2
end

