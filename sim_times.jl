# Julia v1.3
using CALCEPH,PyPlot,Statistics,JLD2,DelimitedFiles,Random,LinearAlgebra,Distributions
rc("font",family="sans-serif")
rc("lines",linewidth=2)
include("regress.jl")
include("CGS.jl")
# Load JPL ephemerides from data and set units
eph = Ephem("../ttv_ss/INPUTS/DE440.bsp") ; prefetch(eph)
options = useNaifId+unitKM+unitDay # useNaifId + unitDay + unitAU
AU = 149597870.700 #km
Random.seed!(42)
# Find when body_id transits between jd1 and jd2 for observer at n_obs with N orbit integration steps 
function find_transit(body_id::Int,eph::CALCEPH.Ephem,jd1::Float64,jd2::Float64,n_obs::Vector{Float64},N::Int)
  JD_0 = 0.0
  ff = Array{Float64}(undef, N)
  xdotn = 0.0
  pos = Array{Float64}(undef,3, N)
  # Compute functions of position and velocity, f(t)=dot(x_bar, v_bar) and f'(t):
  function calc_ffs(t)
    pva = compute(eph,JD_0,t,body_id,10,options,2)./AU
    #println(JD_0)
    x = pva[1:3]; v = pva[4:6]; a = pva[7:9];
    f = dot(x,v) - (dot(x,n_obs))*(dot(v,n_obs))
    Df = dot(v,v) + dot(x,a) - (dot(v,n_obs))^2 - (dot(x,n_obs)*dot(a,n_obs))
    return f,Df,dot(x,n_obs),x
  end
  # Compute minimum sky separation of planet wrt star for all JDs
  dt =  (jd2 - jd1)/(N-1)
  JD = Array{Float64}(undef, N)
  i_min = 1
  ff_min = Inf
  for i=1:N
    JD[i] = jd1 + dt*(i-1)
    JD_0 = JD[i]
    ff[i],Df,xdotn,pos[:,i] = calc_ffs(0.0)
  # Estimate of transit time when f(t)== 0:
      # Df > 0 for transit occuring; 
      # xdotn > 0 for planet in front of star as seen by observer; 
      # local minimum value over entire range (i.e. when close to zero)
    if (Df > 0) && (xdotn > 0) && (abs(ff[i]) < ff_min)
      i_min = i 
      ff_min = abs(ff[i])
    end
  end
   #println("Estimated Transit Time: ",JD[i_min])
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
      f_n,Df_n,xdotn,x = calc_ffs(JD_n)
      JD_n -= f_n/Df_n 
      # Break out if we have reached maximum iterations,or if
      # current transit time estimate equals one of the prior two steps:
      if (JD_n == JD_n1) || (JD_n == JD_n2)
          break
      end
  end          
  JD_tt = JD_0 + JD_n
  #println("Refined Transit Time: ",JD_tt)
  return JD_tt,pos
	# return JD,ff,i_min,pos,JD_tt
end
# Find the transit times for body_id, given planetary period estimate,and number of refinement steps N
function transit_times(body_id::Int,eph::CALCEPH.Ephem,t0,period::Float64,period_err::Float64,n_obs::Vector{Float64},N::Int)
  TT = Float64[]
  nt=1
  # Initialize & find first transit time (with high precision so N=1000):
  JD_tt,pos = find_transit(body_id,eph,t0[1],t0[1]+period,n_obs,1000) 
  push!(TT,JD_tt)
  t_final = t0[end]
  # Find subsequent transit times by shifting time frame by 1 planetary period:
  while JD_tt < t_final # condition to continue shifting frame
    t_start = JD_tt+period-period_err
    t_end = JD_tt+period+period_err
    JD_tt,pos = find_transit(body_id,eph,t_start,t_end,n_obs,N)
    # JD,ff,i_min,pos,JD_tt = find_transit(body_id,eph,t_start,t_end,n_obs,N)
    if (JD_tt>t_final) # last run of while loop doesn't meet condition , so need to break 
      break
    else
      push!(TT,JD_tt)
      nt+=1
    end
  end
  return TT,nt
end
# Add Gaussian sigma noise level (in seconds) to transit times (in days)
function fixed_noise(tt::Vector{Float64},sigma::Real)
  if sigma > 0
    # Draw a random number from a Normal distribution: 
      # dist=Normal(0,sigma / (24 * 3600))
      # noise=rand(dist, length(tt))
    # Since we're assuming that all σ_i are identical:  
      sigtt = ones(length(tt)) * sigma / (24 * 3600) 
      noise = sigtt .* randn(length(tt))
      # println("Noise added with σ of ",string(sigma)," seconds.")
  else
      sigtt=Array{Float64}(undef, length(tt))
      println("No noise added.")
  end
  return sigtt,noise 
end
# Do linear regression of transit times, given mean orbital period
function linear_fit(tt::Vector{Float64},period::Float64,sigtt::Vector{Float64})
  noise = Array{Float64}(undef, length(tt))
  x = Array{Float64}(undef,2, length(tt))
  x[1,1:length(tt)] .= 1.0
  x[2,1] = 0.0 
  for i=2:length(tt)
    # currently accounts for missing transits (noncontinuous) 
    # by rounding [difference in consecutive transit times/Period]
      x[2,i] = round((tt[i]-tt[1])/period) 
  end
  # println(tt,sigtt,std(sigtt))
  # coeff[1] is best linear fit approx of first tt,coeff[2] is average period
  coeff,covcoeff = regress(x,tt,sigtt)
  t0,per=coeff[1],coeff[2]
  return x,t0,per
end
# Collect linear transit times (i.e. t_calc), given mean orbital period
function linear_times(tt::Vector{Float64},period::Float64,sigtt::Vector{Float64})
  x,t0,per=linear_fit(tt,period,sigtt)
  times=collect(t0 .+ per .* range(0,stop = length(tt)-1,length = length(tt))) 
  return times
end
# Compute unit vector which intersects orbital plane of objects 1 and 2:
function calc_obs_loc(pos1,vel1,pos2,vel2)
  h1 = cross(pos1,vel1)
  h2 = cross(pos2,vel2)
  n_obs = cross(h2,h1)
  n_obs /= norm(n_obs) #from one direction when both transit
end
"""
    sim_obs_and_find_times(jd1,sigma,nyear,obs)

 Integrate orbits and find transit times
# Arguments:
- `jd1::Float64`: starting Julian Ephemeris Date of observations.
- `sigma::Real`: fixed noised added to observations.
- `nyear::Real`: time span of observations.
- `obs::String`:source of observations for body 2 (EMB or EV)
# Returns:
   body,tt0,tt,sigtt
"""
function sim_obs_and_find_times(jd1::Float64,sigma::Real,nyear::Real,obs::String)
  # nyear = (jd2 - jd1)/365.25 
  jd2 = nyear*365.25 + jd1
  jdsize = 1000
  # dt = (jd2 - jd1)/jdsize
  @assert (jd1 >= 2287184.5) #2414105.0
  @assert (jd2 <= 2688976.5) #2488985.0
  t0 = range(jd1,stop=jd2,length = jdsize)

  # Compute ephemerides of Sun, Venus and Earth (or EMB):
  pva_sun = Array{Float64}(undef, 9, jdsize)
  pva_venus = Array{Float64}(undef, 9, jdsize)
  pva_earth = Array{Float64}(undef, 9, jdsize)
  for i=1:jdsize
    pva_sun[1:9,i] = compute(eph,t0[i],0.0,10,10,options,2)./AU
    pva_venus[1:9,i] = compute(eph,t0[i],0.0,2,10,options,2)./AU
    if obs=="fromEMB"
      pva_earth[1:9,i] = compute(eph,t0[i],0.0,3,10,options,2)./AU 
    else
      pva_earth[1:9,i] = compute(eph,t0[i],0.0,399,10,options,2)./AU
      # pva_emb = compute(eph,t0[i],0.5,3,10,options,2)
      # pva_moon = compute(eph,t0[i],0.5,301,10,options,2)
      # println("Earth - EMB: ",norm(pva_earth[1:3,i] .- pva_emb[1:3]))
      # println("Earth - Moon: ",norm(pva_earth[1:3,i] .- pva_moon[1:3]))
      # println("Moon - EMB: ",norm(pva_moon[1:3] .- pva_emb[1:3]))
      # println("Ratio: ",norm(pva_earth[1:3,i] .- pva_emb[1:3])/norm(pva_moon[1:3] .- pva_emb[1:3]))
    end
  end

  # Find observer location required to see transits of Venus and Earth:
  n_obs=calc_obs_loc(pva_venus[1:3],pva_venus[4:6],pva_earth[1:3],pva_earth[4:6])

  # Find actual transit times:
  P_venus = 225.0
  P_earth = 365.0
  P_err = 1.0
  tt1,nt1 = transit_times(2,eph,t0,P_venus,P_err,n_obs,10)
  # nt1=length(tt1)
  if obs=="fromEMB"
    tt2,nt2 = transit_times(3,eph,t0,P_earth,P_err,n_obs,10)
  else
    tt2,nt2 = transit_times(399,eph,t0,P_earth,P_err,n_obs,10)
  end
  # nt2=length(tt2)
  # println(nt1," ",nt2)
  # println("Venus Transit Times: ",tt1,'\n',"Earth Transit Times: ",tt2)
	sigtt1,noise1=fixed_noise(tt1,sigma)
	sigtt2,noise2=fixed_noise(tt2,sigma)

  # println("noise to add: ",noise2)
  # for i=1:length(tt2)
  #   println(tt2[i], " ",tt2[i]+noise2[i])
  # end
  tref=2430000
  x1,t01,per1 = linear_fit(tt1+noise1,P_venus,sigtt1)
  x2,t02,per2 = linear_fit(tt2+noise2,P_earth,sigtt2)
  println("P1 linear coefficients: ",t01.-tref," , ",per1)
  println("P2 linear coefficients: ",t02.-tref," , ",per2)
  # println(linear_fit(tt2.+noise2,P_earth,sigtt2))
	tt=[tt1+noise1;tt2+noise2]
  # println("t0= ",t01)
  # println("per= ",per1)
  # # Best-fit linear transit times:
  t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) 
  t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
  # println(t1)
  tt0 = [t1;t2]
	sigtt=[sigtt1;sigtt2]

  body = zeros((nt1+nt2))
  body[1:nt1] .= 1.0
  body[nt1+1:nt1+nt2] .= 2.0

  ttv1=(tt1 .- t1).*24*60
  ttv2=(tt2 .- t2).*24*60
  trans=[round.(range(0,stop = nt1-1, length=nt1));  round.(range(0,stop = nt2-1, length=nt2))]
  function make_transit_times_table()
    name= string("INPUTS/EMBtransit_times",nyear,".txt")
    open(name,"w") do io
      println(io,"# body ntrans  tcalc ttv noise sigma",'\n',"#   JED-2430000 min min min")
      for i=1:nt1
        println(io,"1.0",'\t',i-1,'\t',round(t1[i].-tref,sigdigits=10),'\t',round(ttv1[i],sigdigits=4),'\t',round(noise1[i].*24*60,sigdigits=2),'\t',sigtt1[i].*24*60)
      end
      for i=1:nt2
        println(io,"2.0",'\t',i-1,'\t',round(t2[i].-tref,sigdigits=10),'\t',round(ttv2[i],sigdigits=4),'\t',round(noise2[i].*24*60,sigdigits=2),'\t',sigtt2[i].*24*60)
      end
    end
  end
  make_transit_times_table()
  println("A_TTV1= ",abs((maximum(ttv1)))-abs(minimum(ttv1)))
  println("A_TTV2= ",abs((maximum(ttv2)))-abs(minimum(ttv2)))
  return body,trans,tt,sigtt,tt0
  # return pva_venus,pva_earth
end
# body,tt0,tt,sigtt=sim_obs_and_find_times(2.4332825e6,30,30,"fromEMB")
# Simulate times starting at jd1 for nyear span with sigma Gaussian noise added, save to .txt
function sim_times(jd1::Float64,sigma::Real,nyear::Real,obs::String,dir::String="INPUTS")
  body,tt0,tt,sigtt=sim_obs_and_find_times(jd1,sigma,nyear,obs)
  if obs=="fromEMB"
    name = string(dir,"/EMBtt_",sigma,"s",nyear,"yrs.txt")
  else
    name = string(dir,"/tt_",sigma,"s",nyear,"yrs.txt")
  end
  open(name,"w") do io
    # println(io,"#body",'\t',"tt0",'\t',"tt",'\t',"sigtt")
    for i=1:length(tt)
      println(io,body[i],'\t',tt0[i],'\t',tt[i],'\t',sigtt[i])
    end
  end
end
function sim_times(jd1,nyear,obs)
  body,tt=sim_obs_and_find_times(jd1,0.0,nyear,obs)
  name= string("SS_transit_times.txt")
  open(name,"w") do io
    println(io,"## ",nyear," year long observations starting at ",jd1," JED")
    println(io,"## body_num",'\t',"TT")
  for i=1:length(tt)
    println(io,body[i],'\t',tt[i])
  end
  end  
end

  # Plot orbits along ecliptic and top-down,point to observer of Venus and Earth transits
function plot_orbits(dimension::Int,obs::String,nyear::Real=10)
  jd1=2.4332825e6 ; sigma=30 ;jdsize=1000
  jd2 = nyear*365.25 + jd1
  theta_sun=range(0, stop=2pi, length=100)
  xsun = CGS.RSUN/CGS.AU * cos.(theta_sun)
  ysun = CGS.RSUN/CGS.AU * sin.(theta_sun)
  t0 = range(jd1,stop=jd2-1,length = jdsize)
  pva_sun = zeros(6, jdsize)
  pva_venus = zeros(6, jdsize)
  pva_earth = zeros(6, jdsize)
  pva_emb = zeros(6, jdsize)
  pva_moon = zeros(6, jdsize)
  for i=1:jdsize
    pva_sun[1:6,i] = compute(eph,t0[i],0.0,10,10,options)./AU
    pva_venus[1:6,i] = compute(eph,t0[i],0.0,2,10,options)./AU
    if obs=="fromEMB"
      pva_emb[1:6,i] = compute(eph,t0[i],0.0,3,10,options) ./AU
    else
      pva_moon = compute(eph,t0[i],0.0,301,10,options)./AU
      pva_earth[1:6,i] = compute(eph,t0[i],0.0,399,10,options)./AU
    end
  end
  body,trans,tt,sigtt,tt0=sim_obs_and_find_times(jd1,sigma,nyear,obs)
  nt1 = sum(body .== 1.0)
  nt2 = sum(body .== 2.0)
  trans_pva_venus = zeros(6, jdsize)
  trans_pva_earth = zeros(6, jdsize)
  trans_pva_mars = zeros(6, jdsize)
  trans_pva_jup = zeros(6, jdsize)
  trans_pva_sat = zeros(6, jdsize)
  trans_pva_emb = zeros(6, jdsize)
  trans_pva_moon = zeros(6, jdsize)
  for i=1:length(tt)
    trans_pva_venus[1:6,i] = compute(eph,tt[i],0.0,2,10,options)
    trans_pva_mars[1:6,i] = compute(eph,tt[i],0.0,4,10,options)
    trans_pva_jup[1:6,i] = compute(eph,tt[i],0.0,5,10,options)
    trans_pva_sat[1:6,i] = compute(eph,tt[i],0.0,6,10,options)
    if obs=="fromEMB"
      trans_pva_emb[1:6,i] = compute(eph,tt[i],0.0,3,10,options) 
    else
      trans_pva_moon = compute(eph,tt[i],0.0,301,10,options)
      trans_pva_earth[1:6,i] = compute(eph,tt[i],0.0,399,10,options)
    end

  end

  ## Find position of Moon w.r.t. Earth when Earth transit occurs
  # trans_pva_moon
  

  fig,ax=subplots(figsize=(4,4))#,dpi=150)
  title(string("Location over ",nyear,"yrs"))
  fill(xsun.*5,ysun.*5,color="yellow")
  plot(xsun,ysun,color="yellow")
  plot(pva_venus[1,:],pva_venus[2,:],color="salmon",linewidth=1,alpha=0.5)
  for i=1:nt1
  ax.scatter(trans_pva_venus[1,i],trans_pva_venus[2,i],marker="v",color="salmon",label=string(i))
  end
    if obs=="fromEMB"
      n_obs=calc_obs_loc(trans_pva_venus[1:3],trans_pva_venus[4:6],trans_pva_emb[1:3],trans_pva_emb[4:6])
      plot(pva_emb[1,:],pva_emb[2,:],color="forestgreen",linewidth=1,alpha=0.5)
      ax.scatter(trans_pva_emb[1,nt1+1:nt1+nt2],trans_pva_emb[2,nt1+1:nt1+nt2],marker=".",color="forestgreen",label="EMB")
    else
    n_obs=calc_obs_loc(trans_pva_venus[1:3],trans_pva_venus[4:6],trans_pva_earth[1:3],trans_pva_earth[4:6])
    plot(pva_earth[1,:],pva_earth[2,:],color="forestgreen",linewidth=1,alpha=0.5)
    ax.scatter(trans_pva_earth[1,nt1+1:nt1+nt2],trans_pva_earth[2,nt1+1:nt1+nt2],marker=".",color="forestgreen",label="Earth")
    end
  # arrow(0.0,0.0,n_obs[1],n_obs[2],facecolor="black")
  plot([0,n_obs[1]*1.1],[0,n_obs[2]*1.1],"k--",linewidth=1,alpha=0.5)
  annotate("Line of sight",xy=[n_obs[1];n_obs[2]], xytext=[n_obs[1]+0.05;n_obs[2]],xycoords="data",fontsize="medium") 
  ax.grid(linestyle="--",alpha=0.4)
  xlabel("x [au]",fontsize="large")
  ylabel("y [au]",fontsize="large")
  ylim(-1,1)
  xlim(-1.1,1.1)
  # ax.legend(title="Mid-Transit",fontsize="medium",title_fontsize="medium",markerscale=1.5,loc="upper left")
  # fill(xsun,ysun,color="yellow")
  # plot(xsun,ysun,color="black")
  # plot(pva_venus[1,:],pva_venus[2,:],color="orange",linewidth=1,alpha=0.5)
  # plot(pva_emb[1,:],pva_emb[2,:],color="forestgreen",linewidth=1,alpha=0.5)
  # scatter(trans_pva_venus[1,1:nt1],trans_pva_venus[2,1:nt1],marker=".",color="orange")
  # scatter(trans_pva_emb[1,nt1+1:nt1+nt2],trans_pva_emb[2,nt1+1:nt1+nt2],marker=".",color="forestgreen")
  # plot([0,n_obs[1]*1.2],[0,n_obs[2]*1.2],"k--",linewidth=1)
  # # arrow(0.0,0.0,n_obs[1],n_obs[2],facecolor="black")
  # # ("to observer",xy=[n_obs[1];n_obs[2]], xytext=[n_obs[1]+0.05;n_obs[2]],xycoords="data") 
  # xlabel("x [au]",fontsize=20)
  # ylabel("y [au]",fontsize=20)
  tight_layout()
  @show

  if dimension==3
    fig=figure(figsize=(6,6))
    PyPlot.scatter3D(xsun,ysun,0,marker="o",color=:yellow)
    PyPlot.plot3D(vec(pva_venus[1,:]), vec(pva_venus[2,:]), vec(pva_venus[3,:]),alpha=0.25,color=:orange)
    PyPlot.plot3D(vec(pva_earth[1,:]), vec(pva_earth[2,:]), vec(pva_earth[3,:]),alpha=0.25,color=:forestgreen)
    PyPlot.plot3D([0,n_obs[1]*1.2],[0,n_obs[2]*1.2],[0,n_obs[3]*1.2],linestyle="--",color=:grey)
    # PyPlot.scatter3D(vec(trans_pva_sun[1,:]), vec(trans_pva_sun[2,:]), vec(trans_pva_sun[3,:]),color=:yellow)
    PyPlot.scatter3D(vec(trans_pva_venus[1,1:nt1]),vec(trans_pva_venus[2,1:nt1]),vec(trans_pva_venus[3,1:nt1]),color=:orange,marker=".")
    PyPlot.scatter3D(vec(trans_pva_earth[1,nt1+1:nt1+nt2]),vec(trans_pva_earth[2,nt1+1:nt1+nt2]),vec(trans_pva_earth[3,nt1+1:nt1+nt2]),color=:forestgreen,marker=".")
    PyPlot.tick_params(which="major",
        left=false,right=false,top=false,bottom=false)
    # xlim(-1,1)
    # ylim(-1,1)
    xlabel("x [au]",fontsize=20)
    ylabel("y [au]",fontsize=20)
    zlabel("z [au]",fontsize=20) 
  end
  #   println(n_obs) 
  return trans_pva_venus,trans_pva_emb,trans_pva_mars,trans_pva_jup,trans_pva_sat

  # close()
end
function moon_times(jd1::Float64,sigma::Real,nyear::Real)
  # nyear = (jd2 - jd1)/365.25 
  jd2 = nyear*365.25 + jd1
  jdsize = 1000
  # dt = (jd2 - jd1)/jdsize
  @assert (jd1 >= 2287184.5) #2414105.0
  @assert (jd2 <= 2688976.5) #2488985.0
  t0 = range(jd1,stop=jd2-1,length = jdsize)

  # Compute ephemerides of Sun, Venus and Earth (or EMB):
  pva_sun = zeros(9,jdsize)
  pva_venus = zeros(9,jdsize)
  pva_earth = zeros(9,jdsize)
  #pva_emb = zeros(6, jdsize)
  pva_moon = zeros(9, jdsize)

  for i=1:jdsize
    pva_sun[1:9,i] = compute(eph,t0[i],0.0,10,10,options,2)./AU
    pva_venus[1:9,i] = compute(eph,t0[i],0.0,2,10,options,2)./AU
    pva_earth[1:9,i] = compute(eph,t0[i],0.0,399,10,options,2)./AU
      # pva_emb = compute(eph,t0[i],0.5,3,10,options,2)
    pva_moon[1:9,i] = compute(eph,t0[i],0.5,301,10,options,2)
      # println("Earth - EMB: ",norm(pva_earth[1:3,i] .- pva_emb[1:3]))
      # println("Earth - Moon: ",norm(pva_earth[1:3,i] .- pva_moon[1:3]))
      # println("Moon - EMB: ",norm(pva_moon[1:3] .- pva_emb[1:3]))
      # println("Ratio: ",norm(pva_earth[1:3,i] .- pva_emb[1:3])/norm(pva_moon[1:3] .- pva_emb[1:3]))
    end

  # Find observer location required to see transits of Venus and Earth:
  n_obs=calc_obs_loc(pva_venus[1:3],pva_venus[4:6],pva_earth[1:3],pva_earth[4:6])

  # Find actual transit times:
  P_venus = 225.0
  P_earth = 365.0
  P_err = 1.0
  tt1 = transit_times(2,eph,t0,P_venus,P_err,n_obs,10)
  nt1=length(tt1)
  tt2 = transit_times(399,eph,t0,P_earth,P_err,n_obs,10)
  nt2=length(tt2)
  tt3 = transit_times(301,eph,t0,365.25,P_err,n_obs,10)
  nt3=length(tt3)

  sigtt1=fixed_noise(tt1,sigma)
  sigtt2=fixed_noise(tt2,sigma)
  sigtt3=fixed_noise(tt3,sigma)

  x1,t01,per1 = linear_fit(tt1,P_venus,sigtt1)
  x2,t02,per2 = linear_fit(tt2,P_earth,sigtt2)
  x3,t03,per3 = linear_fit(tt3,365.25,sigtt2)
  println("moon_period=",per3)
  # println("coefficients: ",t02," , ",per2)
  tt=[tt1;tt2;tt3]

  # # Best-fit linear transit times:
  t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) 
  t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
  t3  = collect(t03 .+ per3 .* range(0,stop=nt3-1,length=nt3))
  # # tt0 = [t1;t2;t3]
  # # sigtt=[sigtt1;sigtt2;sigtt3]

  # # body = zeros((nt1+nt2))
  # # body[1:nt1] .= 1.0
  # # body[nt1+1:nt1+nt2] .= 2.0
  subplot(211)
  plot((t2.-t02)./per2,(tt2.-t2).*(24*60))
  xlabel("Time [years]")
  ylabel("TTV [min]")
  subplot(212)
  plot((t3.-t03)./per2,(tt3.-t3).*(24*60))
  xlabel("Time")
  ylabel("TTV [min]")
  toff=((tt2.-t2).*(24*60)).-((tt3.-t3).*(24*60))
  a_s=toff/2pi
  # return tt1,tt2,tt3
  # f=  jldopen(String(FITS/p3moon_fit",sigma,"s",nyear,"yrs.jld2"),"r")
  # tt,tt0,sigtt,ttmodel=f["tt"],f["tt0"],f["sigtt"],f["ttmodel"]
  # pbest_global=f["best_dp"]
  # nplanet,ntrans=f["nplanet"],f["ntrans"]
  # n1,n2=ntrans[1],ntrans[2]
  # mu1,P1,t01,ecos1,esin1=pbest_global[1:5]
  # mu2,P2,t02,ecos2,esin2=pbest_global[6:10]
  # time1=collect(t01 .+ range(0,stop=n1-1,length=n1) .* P1)  #tcalc
  # time2=collect(t02 .+ range(0,stop=n2-1,length=n2) .* P2)
  # tt1,tt2=tt[1:n1],tt[n1+1:n1+n2]       #tobs
  # ttmodel1,ttmodel2 = ttmodel[1:n1],ttmodel[n1+1:n1+n2]
  # ttsim1,ttsim2=(time1.-t01)./365.25,(time2.-t02)./365.25 #in years
  # ttvmodel1,ttvmodel2=(ttmodel1.-time1).*(24*60),(ttmodel2.-time2).*(24*60)
  # ttv1,ttv2=(tt1.-time1).* (24*60),(tt2.-time2).* (24 * 60) #in minutes
  # sigtt1,sigtt2=sigtt[1:n1].* (24 * 60),sigtt[n1+1:n1+n2].* (24 * 60) #in minutes
  # pair_ttvs=decompose_ttvs(f["nplanet"],ntrans[1:f["nplanet"]],pbest_global) .* (24 * 60)

  # moon=moon_ttvs(ntrans,pbest_global) .* (24 * 60)
  # subplot(211)
  # plot(ttsim2,pair_ttvs[2,1,1:n2],color="salmon",label="Venus")
  # plot(ttsim2,moon,linestyle="-.",color="purple",label="Moon")
  # errorbar(ttsim2,ttv2,sigtt2,fmt=".",color="black",capsize=3,ms=5)#,label="Earth")
  # subplot(212)
  # plot(ttsim3,moon,linestyle="-.",color="purple",label="Moon")
  # errorbar(ttsim3,ttv3,sigtt3,fmt=".",color="black",capsize=3,ms=5)#,label="Earth")
  # # text(0,-5.5,label,fontsize="xx-large")
  # xlabel("Time [years]",fontsize=20)
  # ylabel("TTV [min]",fontsize=20)
  # ylim(-7,7)
  # minorticks_on()
  # tick_params(which="both",direction="in",top=true,right=true)
end
# moon_times(jd1,sigma,nyear)

  # function plot_ttvs(sigma)
  #   P_venus = 225
  #   P_earth = 365
  #   P_err = 2
  #   t01,per1=linear_fit(tt1,P_venus,sigma)
  #   t02,per2=linear_fit(tt2,P_venus,sigma)
  #   t1  = linear_times(tt1,P_venus,sigma)
  #   t2  = linear_times(tt2,P_earth,sigma)
  #   # subplot(211)
  #   # scatter((t1.-t01)./per1,tt1.-t1) #x is tranit number 
  #   # plot((t1.- t01)./per1,ttv1) 
  #   # errorbar((t1.-t01)./per1,ttv1,noise1)
  #   scatter((t1.-t01)./365.25,tt1.-t1) # x is JD in years
  #   plot((t1.-t01)./365.25,ttv1)
  #   # subplot(212)
  #   scatter((t2.-t02)./365.25,tt2.-t2,color="green")
  #   plot((t2.-t02)./365.25,ttv2)
  #   errorbar((t2.-t02)./365.25,ttv2,noise2)
  #   # scatter((t2.-t02)./per2,tt2.-t2,color="green") 
  #   # plot((t2.-t02)./per2,ttv2)
  #   # title(sigma)
  #   xlabel("JD (years)")
  #   ylabel("TTVs")
  #   # savefig("OUTPUTs/")
  # end

