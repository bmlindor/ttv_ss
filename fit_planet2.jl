if !@isdefined(TTVFaster)
    include("TTVFaster/src/TTVFaster.jl")
    using Main.TTVFaster
end
import Main.TTVFaster.ttv_wrapper
import Main.TTVFaster.chisquare
include("regress.jl")
using DelimitedFiles,JLD2,Optim,LsqFit,Statistics

function fit_planet2(data_file::Array{}, jmax::Int,tref::Real,tol::Real)
  # (::Core.kwftype(typeof(fit_planet2)))(kws, fit_planet2, data_file, jmax)
  # if haskey(kws, :jmax)
  #     jmax = kws.jmax
  # else
  #     jmax = 5
  # end
  # # etc.
  # Loads .txt datafile with following columns: 
  #body_number (n), initial time of transit (tt0), actual transit times (tt), measurement error (sigtt)
  # data_file = readdlm(filename,Float64)
  nt1 = sum(data_file[:,1] .== 1.0)
  nt2 = sum(data_file[:,1] .== 2.0)
  tt1 = vec(data_file[1:nt1,3]) .- tref
  tt2 = vec(data_file[nt1+1:nt1+nt2,3]) .- tref
  sigtt1 = data_file[1:nt1,4]
  sigtt2 = data_file[nt1+1:nt1+nt2,4]
  # Okay,let's do a linear fit to the transit times (third column):
  function find_coeffs(tt,period,sigtt)
    nt = length(tt)
    x = zeros(2,nt)
    x[1,1:nt] .= 1.0
    x[2,1] = 0.0 
    for i=2:nt
      x[2,i] = x[2,i-1] + round((tt[i]-tt[i-1])/period) 
    end
    coeff,covcoeff = regress(x,tt,sigtt)
    # println(tt,sigtt,std(sigtt))
    println(coeff)
    println(sum((tt .- coeff[1].-coeff[2].*x[2,:])./sigtt).^2)
    return coeff,covcoeff
  end
  # Guess the planets' period by finding median of transit times
  p1_est = median(tt1[2:end] - tt1[1:end-1])
  p2_est = median(tt2[2:end] - tt2[1:end-1])
  # Do fit to find linear ephemeris (t0, per)
  coeff1,covcoeff1 = find_coeffs(tt1,p1_est,sigtt1)
  coeff2,covcoeff2 = find_coeffs(tt2,p2_est,sigtt2)
  t01 = coeff1[1]; per1 = coeff1[2]
  t02 = coeff2[1]; per2 = coeff2[2]
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
    # println("init_param: ",init_param)
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
# If the simulation already exists, can just do 2-planet fit
function fit_planet2(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,obs::String)
  if obs=="fromEMB"
    datafile = string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt")
    outfile = string("FITS/fromEMB/p2_fit",sigma,"s",nyear,"yrs.jld2")
  elseif obs=="fromEV"
    datafile = string("INPUTS/tt_",sigma,"snoEMB",nyear,"yrs.txt")
    outfile = string("FITS/p2_fit",sigma,"s",nyear,"yrs.jld2")
  end
  @assert isfile(datafile)
  println(datafile," loaded.")
  data1 = readdlm(datafile,Float64)
  nt1 = sum(data1[:,1] .== 1.0)
  nt2 = sum(data1[:,1] .== 2.0)
  tt1 = vec(data1[1:nt1,3]) .- tref
  tt2 = vec(data1[nt1+1:nt1+nt2,3]) .- tref
  sigtt1 = data1[1:nt1,4]
  sigtt2 = data1[nt1+1:nt1+nt2,4]

  # Okay,let's do a linear fit to the transit times (third column):
  function find_coeffs(tt,period,sigtt)
    nt = length(tt)
    x = zeros(2,nt)
    x[1,1:nt] .= 1.0
    x[2,1] = 0.0 
    for i=2:nt
      x[2,i] = round((tt[i]-tt[1])/period) 
    end
    coeff,covcoeff = regress(x,tt,sigtt)
    return coeff,covcoeff
  end
  p1est = median(tt1[2:end] - tt1[1:end-1])
  p2est = median(tt2[2:end] - tt2[1:end-1])
  coeff1,covcoeff1 = find_coeffs(tt1,p1est,sigtt1)
  coeff2,covcoeff2 = find_coeffs(tt2,p2est,sigtt2)
  # @assert (sigtt[1] .* (24 * 3600) .= sigma)
  t01 = coeff1[1]; per1 = coeff1[2]
  t02 = coeff2[1]; per2 = coeff2[2]
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
  println("New initial 2-planet fit: ",init_param)

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,init_param)
  best_p2 = fit.param ##### is this the global p2 fit???
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p2,jmax,true)
  lprob_best_p2= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  chi2=chisquare(tt0,nplanet,ntrans,best_p2,tt,sigtt,jmax,true)
  println("Finished 2-planet fit in ",niter," iterations.")
  println("New 2-planet chi-square: ",chi2)
  println("Param: ",best_p2)
  # println("Maximum: ",lprob_best_p2," Param: ",best_p2)
  @save outfile chi2 best_p2 lprob_best_p2 ntrans nplanet tt0 tt ttmodel sigtt
  return best_p2
end
