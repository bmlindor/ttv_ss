if !@isdefined(TTVFaster)
    include("TTVFaster/src/TTVFaster.jl")
    using Main.TTVFaster
end
import Main.TTVFaster.ttv_wrapper
import Main.TTVFaster.chisquare
include("regress.jl")
using DelimitedFiles,JLD2,Optim,LsqFit,Statistics

function fit_moon(filename::String,jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p3in::Float64,p3out::Float64,np3::Int,nphase::Int,dpin::Float64,dpout::Float64,ndp::Int,wide::Bool=false)
  if wide
    fitfile = string("FITS/wide_fit",sigma,"s",nyear,"yrs.jld2")
  else
    fitfile = string("FITS/p3moon_fit",sigma,"s",nyear,"yrs.jld2")
  end
  jd2 = nyear*365.25 + jd1
  data1 = readdlm(filename,Float64)
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
  t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) 
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
  #model = ttv_wrapper2(tt0,param)
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
  #res = optimize(chisquare2,param,method = :l_bfgs,iterations = 21)
  ntrans = [nt1,nt2]
  Nobs = sum(ntrans)
  nplanet = 2
  # create initial simplex? need function for this?
  # result = optimize(f0,xcurr,NelderMead(initial_simplex=MySimplexer(),show_trace=true,iterations=1))
  println("Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,true))
  # res = optimize(params -> chisquare(nplanet,ntrans,params,tt,sigtt),init_param) 
  # init_param = res.minimizer
  param1 = init_param .+ 100.0
  niter = 0
  while maximum(abs.(param1 .- init_param)) > tol && niter < 20
    param1 = init_param
    res = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,init_param)
    init_param = res.param
    niter+=1
    # println("init_param: ",init_param)
    # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt))
  end
  # res = optimize(params -> chisquare(tt0,nplanet,ntrans,params,tt,sigtt),init_param) 
  # init_param = res.minimizer
  # fit2 = curve_fit(ttv_wrapper2,tt0,tt,weight,param; show_trace=true)
  println("Finished 2-planet fit in ",niter," iterations.")
  println("New 2-planet chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,true))
  println("Param: ",init_param)

  # Now,let's add the 3rd planet:
  ntrans = [nt1,nt2,2]
  nplanet = 3
  nparam = 15
  # Grid of periods to search over:
  p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
  p3_cur = 11.86*365.25 
  lprob_p3 = zeros(np3)
  param_p3 = zeros(nparam,np3)
  lprob_best = -1e100 
  p3best = zeros(nparam)
  # Shifting to simulated observation range to search over period grid
  offset = (jd1 + jd2)/2 
  for j=1:np3
    phase = p3[j]*range(0,stop=1,length=nphase) #.+ offset 
    lprob_phase = zeros(nphase)
    lprob_p3[j] = -1e100
    for i=1:nphase 
     # p3 param_names: mass ratio,phase,ecosw,esinw
      param_tmp = [log10(1e-3),phase[i],0.01,0.01] 
      param3 = [init_param;param_tmp] 
      p3_cur = p3[j] 
      # fit = curve_fit(ttv_wrapper_fixp3,tt0,tt,weight,param3) #optimizes fit w/ 3 planet model
      # fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,true,p3_cur),tt0,tt,weight,param3) 
      # param3 = fit.param
      param1 = param3 .+ 100.0
      niter = 0
      while maximum(abs.(param1 .- param3)) > tol && niter < 20
        param1 = param3
        fit = curve_fit((tt0,param3) -> ttv_wrapper(tt0,nplanet,ntrans,[param3[1:10];10^param3[11];p3_cur;param3[12:end]],jmax,true),tt0,tt,weight,param3)
        param3 = fit.param
        niter+=1
        # println("init_param: ",param3)
        # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,param3,tt,sigtt,true,p3_cur))
      end
      ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[param3[1:10];10^param3[11];p3_cur;param3[12:end]],jmax,true)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best 
        lprob_best = lprob_phase[i]
        p3best = [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      end
      if lprob_phase[i] > lprob_p3[j] 
        lprob_p3[j] = lprob_phase[i]
        param_p3[1:nparam,j] =  [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      end
    end
    # println("Period: ",p3[j]," log Prob: ",lprob_p3[j]," Param: ",vec(param_p3[1:nparam,j]))
  end
  println("Finished 3-planet fit w/ fixed period: ",p3best," in ",niter," iterations")

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,p3best)
  best_p3 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p3,jmax,true)
  lprob_best_p3= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  println("Finished global 3-planet fit.")
  println("New 3-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p3,tt,sigtt,jmax,true))
  println("Maximum: ",lprob_best_p3," Param: ",best_p3)

  # Now,search for Moon:
  nparam = 18
  deltaphi_cur = 2.312
  deltaphi = range(dpin,stop=dpout,length=ndp)
  lprob_dp = zeros(ndp)
  param_dp = zeros(nparam,ndp)
  lprob_best = -1e100 
  dpbest = zeros(nparam)
  niter = 0
  for j=1:ndp
    lprob_dp[j] = -1e100 
    # lunar params: t_s ,t_c ,deltaphi 
    param_tmp = [0.01,0.01,deltaphi[j]] 
    param5 = [best_p3;param_tmp]
    deltaphi_cur = deltaphi[j]
    param1 = param5 .+ 100.0
    while maximum(abs.(param1 .- param5)) > tol && niter < 20
      param1 = param5
      fit = curve_fit((tt0,param5) -> ttv_wrapper(tt0,nplanet,ntrans,param5,jmax,false),tt0,tt,weight,param5)
      param5 = fit.param 
      niter+=1
    end
    ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[fit.param[1:17];deltaphi_cur],jmax,false)
    lprob_dp[j]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
    if lprob_dp[j] > lprob_best 
      lprob_best = lprob_dp[j]
      dpbest = [fit.param[1:17];deltaphi_cur]
    end
    param_dp[1:nparam,j] = [fit.param[1:17];deltaphi_cur]
    # println("deltaphi: ",deltaphi[j]," log Prob: ",lprob_dp[j]," Param: ",vec(param_dp[1:nparam,j]))
  end
  println("Finished lunar search: ",dpbest," in ",niter," iterations")

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,false),tt0,tt,weight,dpbest)
  best_dp = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_dp,jmax,false)
  lprob_best_dp = (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  println("Finished global moon fit.")
  println("3-planet lunar chi-square: ",chisquare(tt0,nplanet,ntrans,best_dp,tt,sigtt,jmax,false))
  println("Maximum: ",lprob_best_dp," Param: ",best_dp)

  @save fitfile p3 lprob_p3 best_p3 lprob_best_p3 deltaphi lprob_dp best_dp lprob_best_dp ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase dpin dpout ndp
  return best_p3,best_dp
end
# If planet fit already exists, can just do moon search
function fit_moon(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,dpin::Float64,dpout::Float64,ndp::Int,nplanets::Real)
  infile = string("FITS/p",nplanets,"_fit",sigma,"s",nyear,"yrs.jld2")
  outfile = string("FITS/p",nplanets,"moon_fit",sigma,"s",nyear,"yrs.jld2")
  @assert isfile(infile)
  m = jldopen(String(infile),"r")
  tt0,tt,ttmodel,sigtt=m["tt0"],m["tt"],m["ttmodel"],m["sigtt"]
  nt1,nt2 = m["ntrans"][1],m["ntrans"][2]
  nplanet,ntrans = m["nplanet"],m["ntrans"]
  if nplanets == 2
    best_par = m["init_param"]
  else 
    per, lprob_per = m[string("p",nplanets)],m[string("lprob_p",nplanets)]
    best_par = m[string("best_p",nplanets)]
    lprob_best_par = m[string("lprob_best_p",nplanets)]
  end
  Nobs = sum([nt1,nt2])
  jmax=5
  jd2 = nyear*365.25 + jd1
  weight = ones(nt1+nt2)./ sigtt.^2
  println(infile," loaded.")
  println("Previous model params: ",best_par)

  # Now,search for Moon:
  nparam = length(best_par)+3
  deltaphi_cur = 2.312
  deltaphi = range(dpin,stop=dpout,length=ndp)
  lprob_dp = zeros(ndp)
  param_dp = zeros(nparam,ndp)
  lprob_best = -1e100 
  dpbest = zeros(nparam)
  niter = 0
  for j=1:ndp
    lprob_dp[j] = -1e100 
    # lunar params: tmax sin(phi_0) ,tmax cos(phi_0) ,deltaphi 
    param_tmp = [0.01,0.01,deltaphi[j]] 
    param5 = [best_par;param_tmp]
    deltaphi_cur = deltaphi[j]
    param1 = param5 .+ 100.0
    while maximum(abs.(param1 .- param5)) > tol #&& niter < 20
      param1 = param5
      fit = curve_fit((tt0,param5) -> ttv_wrapper(tt0,nplanet,ntrans,param5,jmax,false),tt0,tt,weight,param5)
      param5 = fit.param
      niter+=1 
    end
    ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[fit.param[1:end-1];deltaphi_cur],jmax,false)
    lprob_dp[j]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
    if lprob_dp[j] > lprob_best 
      lprob_best = lprob_dp[j]
      dpbest = [fit.param[1:end-1];deltaphi_cur]
    end
    param_dp[1:nparam,j] = [fit.param[1:end-1];deltaphi_cur]
    # println("deltaphi: ",deltaphi[j]," log Prob: ",lprob_dp[j]," Param: ",vec(param_dp[1:nparam,j]))
  end
  println("Finished lunar search: ",dpbest," in ",niter," iterations")

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,false),tt0,tt,weight,dpbest)
  best_dp = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_dp,jmax,false)
  lprob_best_dp = (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  println("Finished global moon fit.")
  println("Lunar chi-square: ",chisquare(tt0,nplanet,ntrans,best_dp,tt,sigtt,jmax,false))
  println("Maximum: ",lprob_best_dp," Param: ",best_dp)
  @save outfile deltaphi lprob_dp best_dp lprob_best_dp ntrans nplanet tt0 tt ttmodel sigtt
  return best_dp
end
