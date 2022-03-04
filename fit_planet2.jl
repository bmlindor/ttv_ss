if !@isdefined(TTVFaster)
    # include("TTVFaster/src/TTVFaster.jl")
  # using TTVFaster
end
import TTVFaster.ttv_wrapper
import TTVFaster.chisquare
include("regress.jl")
using DelimitedFiles,JLD2,Optim,LsqFit,Statistics
# If the simulation already exists, can just do 2-planet fit
function fit_planet2(
  jd1::Float64,sigma::Float64,nyear::Float64,
  addnoise::Bool=true,EM::Bool=true)
  if EM
    datafile = string("INPUTS/tt_",sigma,"sEMB",nyear,"yrs.txt")
  else
    datafile = string("INPUTS/tt_",sigma,"snoEMB",nyear,"yrs.txt")
  end
  jd2 = nyear*365.25 + jd1
  data1 = readdlm(datafile,Float64)
  nt1 = sum(data1[:,1] .== 1.0)
  nt2 = sum(data1[:,1] .== 2.0)
  tt1 = vec(data1[1:nt1,3])
  tt2 = vec(data1[nt1+1:nt1+nt2,3])
  if addnoise 
    sigtt1 = data1[1:nt1,4]
    sigtt2 = data1[nt1+1:nt1+nt2,4]
  else
    sigtt1 = ones(nt1)
    sigtt2 = ones(nt2)
  end
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
    return coeff,covcoeff
  end
  p1est = median(tt1[2:end] - tt1[1:end-1])
  p2est = median(tt2[2:end] - tt2[1:end-1])
  coeff1,covcoeff1 = find_coeffs(tt1,p1est,sigtt1)
  coeff2,covcoeff2 = find_coeffs(tt2,p2est,sigtt2)
  sigtt=[sigtt1;sigtt2] 
  # @assert (sigtt[1] .* (24 * 3600) .= sigma)
  t01 = coeff1[1]; per1 = coeff1[2]
  t02 = coeff2[1]; per2 = coeff2[2]
  t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) #best fit linear transit times w/o ttvs
  t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
  # Best-fit linear transit times:
  tt0 = [t1;t2]
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
  # Actual transit times:
  tt=[tt1;tt2]
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
  println("Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,EM))
  param1 = init_param .+ 100.0
  while maximum(abs.(param1 .- init_param)) > 1e-5
    param1 = init_param
    res = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,init_param)
    init_param = res.param
    # println("init_param: ",init_param)
    # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt))
  end
  println("Finished 2-planet fit: ",init_param)
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,init_param,jmax,EM)
  lprob_best_p2= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  chi2=chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,EM)
  println("New 2-planet chi-square: ",chi2)
  if EM
    fitfile = string("FITS/fromEMB/p2_fit",sigma,"s",nyear,"yrs.jld2")
  else
    fitfile = string("FITS/p2_fit",sigma,"s",nyear,"yrs.jld2")
  end
  @save fitfile chi2 init_param ntrans nplanet tt0 tt ttmodel sigtt
  return 
end
