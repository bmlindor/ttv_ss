include("sim_times.jl")
include("CGS.jl")
using DelimitedFiles,JLD2,Optim,LsqFit,Statistics

function fit_planet2(filename::String, jmax::Int,tref::Real,tol::Real)
  # (::Core.kwftype(typeof(fit_planet2)))(kws, fit_planet2, data_file, jmax)
  # if haskey(kws, :jmax)
  #     jmax = kws.jmax
  # else
  #     jmax = 5
  # end
  # # etc.
  # Loads .txt datafile with following columns: 
  #body_number (n), initial time of transit (tt0), actual transit times (tt), measurement error (sigtt)
  data_file = readdlm(filename,Float64)
  nt1 = sum(data_file[:,1] .== 1.0)
  nt2 = sum(data_file[:,1] .== 2.0)
  tt1 = vec(data_file[1:nt1,3]) .- tref
  tt2 = vec(data_file[nt1+1:nt1+nt2,3]) .- tref
  sigtt1 = data_file[1:nt1,4]
  sigtt2 = data_file[nt1+1:nt1+nt2,4]
  # Okay,let's do a linear fit to the transit times (third column):
  # Guess the planets' period by finding median of transit times
  p1est = median(tt1[2:end] - tt1[1:end-1])
  p2est = median(tt2[2:end] - tt2[1:end-1])
  x1,t01,per1 = linear_fit(tt1,p1est,sigtt1)
  x2,t02,per2 = linear_fit(tt2,p2est,sigtt2)
  t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) 
  t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
  # Best-fit linear transit times:
  tt0 = [t1;t2]
  # Actual transit times:
  tt=[tt1;tt2]
  sigtt=[sigtt1;sigtt2] 
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2

  # Okay,now let's do a 2-planet fit:
  # param_names = mass ratio,period,initial transit time,e*cos(omega),e*sin(omega)
  init_param = [3e-6,per1,t01,0.01,0.01,
                3e-6,per2,t02,0.01,0.01] 
  println("Initial parameters: ",init_param)
  # Set up data structure to hold planet properties,passed to TTVFaster
  data=init_param
  p1=TTVFaster.Planet_plane_hk(data[1],data[2],data[3],data[4],data[ 5])
  p2=TTVFaster.Planet_plane_hk(data[6],data[7],data[8],data[9],data[10])
  time1 = collect(p1.trans0 .+ range(0,stop=nt1-1,length=nt1) .* p1.period)
  time2 = collect(p2.trans0 .+ range(0,stop=nt2-1,length=nt2) .* p2.period)
  # Initialize the computation of the Laplace coefficients:
  ttv1 = zeros(nt1)
  ttv2 = zeros(nt2)
  # Need first call to TTVFaster,without optimizing
  dummy=TTVFaster.compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2)

  # Now,optimize 2-planet fit
  ntrans = [nt1,nt2]
  Nobs = sum(ntrans)
  nplanet = 2
  println("Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,true))
  param1 = init_param .+ 100.0
  niter = 0
  # pscale = [1e-5,1.0,1.0,0.01,0.01,1e-5,1.0,1.0,0.01,0.01]
  pscale = ones(10)
  while maximum(abs.(param1 .- init_param)./pscale) > tol
    param1 = init_param
    fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,init_param)
    init_param = fit.param
    # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,true))
    niter += 1
  end
  println("Finished 2-planet fit: ",init_param," ",niter)
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,init_param,jmax,true)
  lprob_best= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  chi2=chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,true)
  println("New 2-planet chi-square: ",chi2) 
  println("Output file in Fitfromdatafile.jld2")
  outfile = string("Fitfromdatafile.jld2")
  @save outfile chi2 init_param lprob_best ntrans nplanet tt0 tt ttmodel sigtt
  return 
end

"""
    fit_planet2(jd1,sigma,nyear,tref,tol,obs)

 If the observations already exists, do 2-planet fit.
# Arguments:
- `jd1::Float64`: starting Julian Ephemeris Date of observations.
- `sigma::Real`: fixed noised added to observations.
- `nyear::Real`: time span of observations.
- `tref::Real`: JED to subtract from transit times to aid fit of low mass planets.
- `tol::Real`: tolerance level of fit.
- `obs::String`: source of observations for body 2 (EMB or EV).
# Returns:
- `best_p2::Vector{Float64}`: list of global best paramters for 2 planets given the observed transit times.
- `lprob_best_p2::Float64`: log probability of detecting 2 planets with the given properties.
"""
function fit_planet2(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,obs::String)
  if obs=="fromEMB"
    datafile = string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt")
    outfile = string("FITS/fromEMB/p2_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/fromEMB/p2_fit",sigma,"s",nyear,"yrs.txt")
  elseif obs=="fromEV"
    datafile = string("INPUTS/tt_",sigma,"snoEMB",nyear,"yrs.txt")
    outfile = string("FITS/p2_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/p2_fit",sigma,"s",nyear,"yrs.txt")
  end
  @assert isfile(datafile)
  println(datafile," loaded.")
  data1 = readdlm(datafile,Float64,comments=true)
  nt1 = sum(data1[:,1] .== 1.0)
  nt2 = sum(data1[:,1] .== 2.0)
  tt1 = vec(data1[1:nt1,3]) .- tref
  tt2 = vec(data1[nt1+1:nt1+nt2,3]) .- tref
  sigtt1 = data1[1:nt1,4]
  sigtt2 = data1[nt1+1:nt1+nt2,4]

  # Okay,let's do a linear fit to the transit times (third column):
  p1est = median(tt1[2:end] - tt1[1:end-1])
  p2est = median(tt2[2:end] - tt2[1:end-1])
  x1,t01,per1 = linear_fit(tt1,p1est,sigtt1)
  x2,t02,per2 = linear_fit(tt2,p2est,sigtt2)
  t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) 
  t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
  # Best-fit linear transit times without TTVs:
  tt0 = [t1;t2]
  # Actual transit times:
  tt=[tt1;tt2]
  sigtt=[sigtt1;sigtt2]
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
  # Okay,now let's do a 2-planet fit:
  # param_names = mass ratio,period,initial transit time,e*cos(omega),e*sin(omega)
  init_param = [3e-6,per1,t01,0.01,0.01,
                3e-6,per2,t02,0.01,0.01] 
  println("Initial parameters: ",init_param)
  # Set up data structure to hold planet properties,passed to TTVFaster
  jmax = 5
  data=init_param
  p1=TTVFaster.Planet_plane_hk(data[1],data[2],data[3],data[4],data[ 5])
  p2=TTVFaster.Planet_plane_hk(data[6],data[7],data[8],data[9],data[10])
  # assuming no transits are skipped/duplicated ########### To Change
  time1 = collect(p1.trans0 .+ range(0,stop=nt1-1,length=nt1) .* p1.period)
  time2 = collect(p2.trans0 .+ range(0,stop=nt2-1,length=nt2) .* p2.period)
  ##########
  # Initialize the computation of the Laplace coefficients:
  ttv1 = zeros(nt1)
  ttv2 = zeros(nt2)
  # Need first call to TTVFaster,without optimizing
  dummy=TTVFaster.compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2)

  # Now,optimize 2-planet fit
  ntrans = [nt1,nt2]
  Nobs = sum(ntrans)
  nplanet = 2
	nparam=10
  println("Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,true))
  param1 = init_param .+ 100.0
  niter = 0
  while maximum(abs.(param1 .- init_param)) > tol && niter < 20
    param1 = init_param
    res = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,init_param)
    init_param = res.param
    niter += 1
    # println("init_param: ",init_param)
    # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt))
  end
  println("New initial 2-planet fit: ",init_param," in ",niter," iterations.")

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,init_param)
  cov=estimate_covar(fit)
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  best_p2 = fit.param ##### is this the global p2 fit???
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p2,jmax,true)
  lprob_best_p2= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished 2-planet fit") 
  println("New 2-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p2,tt,sigtt,jmax,true))
  println("Maximum: ",lprob_best_p2," Param: ",best_p2)

  pname=["mu_1","P_1","t01","ecos1","esin1",
          "mu_2","P_2","t02","ecos2","esin2"]
  mean_mp=[best_p2[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_p2[(iplanet-1)*5+4]^2 + best_p2[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  ecc_errs=[sqrt(err[(iplanet-1)*5+4]^2 + err[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]

  open(results,"w") do io
    println(io,"Global Fit Results.")
    for i=1:nparam
      println(io,pname[i],": ",best_p2[i]," ± ",err[i])
    end
    println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",'\n',mean_ecc,'\n'," ± ",ecc_errs)
  end
   @save outfile best_p2 lprob_best_p2 ntrans nplanet tt0 tt ttmodel sigtt
  return best_p2,lprob_best_p2
end
