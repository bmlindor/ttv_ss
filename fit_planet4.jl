if !@isdefined(TTVFaster)
    include("TTVFaster/src/TTVFaster.jl")
    using Main.TTVFaster
end
import Main.TTVFaster.ttv_wrapper
import Main.TTVFaster.chisquare
include("regress.jl")
using DelimitedFiles,JLD2,Optim,LsqFit,Statistics

function fit_planet4(filename::String,jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p3in::Float64,p3out::Float64,np3::Int,nphase::Int,p4in::Float64,p4out::Float64,np4::Int,obs::String)
  if obs=="fromEMB"
    fitfile = string("FITS/fromEMB/p4_fit",sigma,"s",nyear,"yrs.jld2")
  elseif obs=="fromEV"
    fitfile = string("FITS/p4_fit",sigma,"s",nyear,"yrs.jld2")
  end
  jd2 = nyear*365.25 + jd1
  data1 = readdlm(filename)
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
  println("Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,EM))
  # res = optimize(params -> chisquare(nplanet,ntrans,params,tt,sigtt),init_param) 
  # init_param = res.minimizer
  param1 = init_param .+ 100.0
  niter = 0
  while maximum(abs.(param1 .- init_param)) > tol && niter < 20
    param1 = init_param
    res = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,init_param)
    init_param = res.param
    niter+=1
    # println("init_param: ",init_param)
    # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt))
  end
  # res = optimize(params -> chisquare(tt0,nplanet,ntrans,params,tt,sigtt),init_param) 
  # init_param = res.minimizer
  # fit2 = curve_fit(ttv_wrapper2,tt0,tt,weight,param; show_trace=true)
  println("Finished 2-planet fit in ",niter," iterations.")
  println("New p2 chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,EM))
  println("Param: ",init_param)

  # Now,let's add the 3rd planet:
  ntrans = [nt1,nt2,2] 
  nplanet = 3
  nparam = 15
  # Grid of periods to search over:
  p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
  lprob_p3 = zeros(np3)
  p3_cur = 11.86*365.25 #jupiter period in days,initial value
  param_p3 = zeros(nparam,np3)
  lprob_best = -1e100 #global best fit
  p3best = zeros(nparam)
  niter = 0
  for j=1:np3
    phase = p3[j]*range(0,stop=1,length=nphase) 
    lprob_phase = zeros(nphase)
    lprob_p3[j] = -1e100
    for i=1:nphase #loops over jupiter phases
     # p3 param_names: mass ratio,phase,ecosw,esinw
      param_tmp = [log10(1e-3),phase[i],0.01,0.01] 
      param3 = [init_param;param_tmp] #concatenate 2 planet model to 3 planet model params
      p3_cur = p3[j] #sets jupiter period to global value
      param1 = param3 .+ 100.0
      niter=0
      while maximum(abs.(param1 .- param3)) > tol && niter < 20
        param1 = param3
        fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:10];10^params[11];p3_cur;params[12:end]],jmax,EM),tt0,tt,weight,param3)
        param3 = fit.param
        niter+=1
        # println("init_param: ",param3)
        # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,param3,tt,sigtt,true,p3_cur))
      end
      ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[param3[1:10];10^param3[11];p3_cur;param3[12:end]],jmax,EM)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best # check that best fit for period is better than global best fit
        lprob_best = lprob_phase[i]
        p3best = [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      end
      if lprob_phase[i] > lprob_p3[j] # checks best fit over all phases of jupiter for this particular period
        lprob_p3[j] = lprob_phase[i]
        param_p3[1:nparam,j] =  [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      end
    end
    # println("Period: ",p3[j]," log Prob: ",lprob_p3[j]," Param: ",vec(param_p3[1:nparam,j]))
  end
  println("Finished 3-planet fit w/ fixed period: ",p3best," in ",niter," iterations")
  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,p3best)
  best_p3 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p3,jmax,EM)
  lprob_best_p3= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  println("Finished global 3-planet fit.")
  println("New 3-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p3,tt,sigtt,jmax,EM))
  println("Maximum: ",lprob_best_p3," Param: ",best_p3)

  # Now,add a 4th planet:
  ntrans = [nt1,nt2,2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
  nplanet = 4
  nparam = 20
  # Grid of periods to search over:
  p4 = 10 .^ range(log10(p4in),stop=log10(p4out),length=np4)
  p4_cur =  1.88*365.25 
  lprob_p4 = zeros(np4)
  param_p4 = zeros(nparam,np4)
  lprob_best = -1e100 
  p4best = zeros(nparam)
  niter = 0
  for j=1:np4
    phase = p4[j]*range(0,stop=1,length=nphase)  
    lprob_phase = zeros(nphase) 
    lprob_p4[j] = -1e100
    for i=1:nphase
     # p4 param_names: mass ratio,phase,ecosw,esinw; uses same nphase as p3
      param_tmp = [log10(1e-7),phase[i],0.01,0.01]
      # Mars' period is shorter than Jupiter's, so need to keep sorted for now
      param4 = [best_p3[1:10];param_tmp;best_p3[11:15]]   
      p4_cur = p4[j]
      param1 = param4 .+ 100.0
      niter=0
      while maximum(abs.(param1 .- param4)) > tol #&& niter < 20
        param1 = param4
        fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:10];10^params[11];p4_cur;params[12:end]],jmax,EM),tt0,tt,weight,param4)
        param4 = fit.param 
        niter+=1
        # println("init_param: ",param4)
        # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,param4,tt,sigtt,5,true))
      end
      ttmodel=ttv_wrapper(tt0,nplanet,ntrans,[param4[1:10];10^param4[11];p4_cur;param4[12:end]],jmax,EM)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best
        lprob_best = lprob_phase[i]
        p4best = [fit.param[1:10];10^param4[11];p4_cur;fit.param[12:end]]
      end
      if lprob_phase[i] > lprob_p4[j] 
        lprob_p4[j] = lprob_phase[i]
        param_p4[1:nparam,j] =  [fit.param[1:10];10^param4[11];p4_cur;fit.param[12:end]]
      end
    end
    # println("Period: ",p4[j]," log Prob: ",lprob_p4[j]," Param: ",vec(param_p4[1:nparam,j]))
  end
  println("Finished 3-planet fit w/ fixed period: ",p3best," in ",niter," iterations")

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,p4best)
  best_p4 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p4,jmax,EM)
  lprob_best_p4= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  println("Finished global 4-planet fit.")
  println("New 4-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p4,tt,sigtt,jmax,EM))
  println("Maximum: ",lprob_best_p4," Param: ",best_p4)

  @save fitfile p3 lprob_p3 best_p3 lprob_best_p3 p4 lprob_p4 best_p4 lprob_best_p4 ntrans nplanet tt0 tt ttmodel sigtt
  return best_p3,best_p4 
end
# If 3-planet fit already exists, can just do 4-planet search
function fit_planet4(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p4in::Float64,p4out::Float64,np4::Int,nphase::Int,obs::String)
  if obs=="fromEMB"
    infile = string("FITS/fromEMB/p3_fit",sigma,"s",nyear,"yrs.jld2")
    outfile = string("FITS/fromEMB/p4_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/fromEMB/p4_fit",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/fromEMB/p4_grid",sigma,"s",nyear,"yrs.txt")
  elseif obs=="fromEV"
    infile = string("FITS/p3_fit",sigma,"s",nyear,"yrs.jld2")
    outfile = string("FITS/p4_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/p4_fit",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/p4_grid",sigma,"s",nyear,"yrs.txt")
  end
  @assert isfile(infile)
  m = jldopen(String(infile),"r")
  tt0,tt,ttmodel,sigtt=m["tt0"],m["tt"],m["ttmodel"],m["sigtt"]
  nt1,nt2 = m["ntrans"][1],m["ntrans"][2]
  p3,lprob_p3=m["p3"],m["lprob_p3"]
  best_p3,lprob_best_p3=m["best_p3"],m["lprob_best_p3"]
  Nobs = sum([nt1,nt2])
  jmax=5
  jd2 = nyear*365.25 + jd1
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
  println(infile," loaded.")
  println("Previous model params: ",best_p3)

 # Now,add a 4th planet:
  ntrans = [nt1,nt2,2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
  nplanet = 4
  nparam = 20
  # Grid of periods to search over:
  p4 = 10 .^ range(log10(p4in),stop=log10(p4out),length=np4)
  p4_cur =  1.88*365.25 
  lprob_p4 = zeros(np4)
  param_p4 = zeros(nparam,np4)
  lprob_best = -1e100 
  p4best = zeros(nparam)
  niter = 0
  for j=1:np4
    phase = p4[j]*range(0,stop=1,length=nphase)  
    lprob_phase = zeros(nphase) 
    lprob_p4[j] = -1e100
    for i=1:nphase
     # p4 param_names: mass ratio,phase,ecosw,esinw; uses same nphase as p3
      param_tmp = [log10(1e-7),phase[i],0.01,0.01]
      # Mars' period is shorter than Jupiter's, so need to keep sorted for now
      param4 = [best_p3[1:10];param_tmp;best_p3[11:15]]   
      p4_cur = p4[j]
      param1 = param4 .+ 100.0
      niter=0
      while maximum(abs.(param1 .- param4)) > tol && niter < 20
        param1 = param4
        fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:10];10^params[11];p4_cur;params[12:end]],jmax,true),tt0,tt,weight,param4)
        param4 = fit.param 
        niter+=1
      end
      ttmodel=ttv_wrapper(tt0,nplanet,ntrans,[param4[1:10];10^param4[11];p4_cur;param4[12:end]],jmax,true)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best
        lprob_best = lprob_phase[i]
        p4best = [fit.param[1:10];10^param4[11];p4_cur;fit.param[12:end]]
      end
      if lprob_phase[i] > lprob_p4[j] 
        lprob_p4[j] = lprob_phase[i]
        param_p4[1:nparam,j] =  [fit.param[1:10];10^param4[11];p4_cur;fit.param[12:end]]
      end
    end
    # println("Period: ",p4[j]," log Prob: ",lprob_p4[j]," Param: ",vec(param_p4[1:nparam,j]))
  end
  println("Finished 4-planet fit w/ fixed period: ",p4best," in ",niter," iterations")
  writedlm(grid,zip(p4,lprob_p4))
  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,p4best)
  cov=estimate_covar(fit)
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  best_p4 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p4,jmax,true)
  lprob_best_p4= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global 4-planet fit.")
  println("New 4-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p4,tt,sigtt,jmax,true))
  println("Maximum: ",lprob_best_p4," Param: ",best_p4)

  pname=["mu_1","P_1","t01","ecos1","esin1",
        "mu_2","P_2","t02","ecos2","esin2",
        "mu_3","P_3","t03","ecos3","esin3",
        "mu_4","P_4","t04","ecos4","esin4"]
  mean_mp=[best_p4[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_p4[(iplanet-1)*5+4]^2 + best_p4[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  ecc_errs=[sqrt(err[(iplanet-1)*5+4]^2 + err[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]

  open(results,"w") do io
    println(io,"Global Fit Results. Per=[",p4in," - ",p4out,", length=",np4,"]")
    for i=1:nparam
      println(io,pname[i],": ",best_p4[i]," ± ",err[i])
    end
    println(io,"Retrieved Earth masses: ",mean_mp," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",mean_ecc," ± ",ecc_errs)
  end
  @save outfile p3 lprob_p3 best_p3 lprob_best_p3 p4 lprob_p4 best_p4 lprob_best_p4 ntrans nplanet tt0 tt ttmodel sigtt nphase
  return best_p4 
end
# If 3-planet fit with moon already exists, can do 4-planet search
function fit_planet4(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p4in::Float64,p4out::Float64,np4::Int,nphase::Int)
  infile = string("FITS/p3moon_fit",sigma,"s",nyear,"yrs.jld2")
  outfile = string("FITS/moonp4_fit",sigma,"s",nyear,"yrs.jld2")
  @assert isfile(infile)
  m = jldopen(String(infile),"r")
  tt0,tt,ttmodel,sigtt=m["tt0"],m["tt"],m["ttmodel"],m["sigtt"]
  nt1,nt2 = m["ntrans"][1],m["ntrans"][2]
  p3,lprob_p3=m["p3"],m["lprob_p3"]
  best_p3,lprob_best_p3=m["best_p3"],m["lprob_best_p3"]
  dpin,dpout,ndp = m["dpin"],m["dpout"],m["ndp"]
  dp =  range(dpin,stop=dpout,length=ndp)
  lprob_dp=m["lprob_dp"]
  best_dp,lprob_best_dp=m["best_dp"],m["lprob_best_dp"]
  Nobs = sum([nt1,nt2])
  jmax=5
  jd2 = nyear*365.25 + jd1
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
  println(infile," loaded.")
  println("Previous model params: ",best_dp)
 # Now,add a 4th planet:
  ntrans = [nt1,nt2,2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
  nplanet = 4
  nparam = 23
  # Grid of periods to search over:
  p4 = 10 .^ range(log10(p4in),stop=log10(p4out),length=np4)
  p4_cur =  1.88*365.25 
  lprob_p4 = zeros(np4)
  param_p4 = zeros(nparam,np4)
  lprob_best = -1e100 
  p4best = zeros(nparam)
  niter = 0
  for j=1:np4
    phase = p4[j]*range(0,stop=1,length=nphase) 
    lprob_phase = zeros(nphase) 
    lprob_p4[j] = -1e100
    for i=1:nphase
     # p4 param_names: mass ratio,phase,ecosw,esinw; uses same nphase as p3
      param_tmp = [log10(1e-7),phase[i],0.01,0.01]
      # Mars' period is shorter than Jupiter's, so need to keep sorted for now
      param4 = [best_dp[1:10];param_tmp;best_dp[11:end]]   
      p4_cur = p4[j]
      param1 = param4 .+ 100.0
      niter=0
      while maximum(abs.(param1 .- param4)) > tol && niter <20
        param1 = param4
        # println("param1 ", param1)
        fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:10];10^params[11];p4_cur;params[12:end]],jmax,false),tt0,tt,weight,param4)
        param4 = fit.param 
        niter+=1
        # println(p4_cur)
        println("param4 ", param4)
      end
      ttmodel=ttv_wrapper(tt0,nplanet,ntrans,[param4[1:10];10^param4[11];p4_cur;param4[12:end]],jmax,false)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best
        lprob_best = lprob_phase[i]
        p4best = [fit.param[1:10];10^param4[11];p4_cur;fit.param[12:end]]
      end
      if lprob_phase[i] > lprob_p4[j] 
        lprob_p4[j] = lprob_phase[i]
        param_p4[1:nparam,j] =  [fit.param[1:10];10^param4[11];p4_cur;fit.param[12:end]]
      end
    end
    # println("Period: ",p4[j]," log Prob: ",lprob_p4[j]," Param: ",vec(param_p4[1:nparam,j]))
  end
  println("Finished 4-planet fit w/ fixed period: ",p4best," in ",niter," iterations")
  writedlm(grid,zip(p4,lprob_p4))
  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,false),tt0,tt,weight,p4best)
  cov=estimate_covar(fit)
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  best_p4 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p4,jmax,false)
  lprob_best_p4= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global 4-planet fit.")
  println("New 4-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p4,tt,sigtt,jmax,false))
  println("Maximum: ",lprob_best_p4," Param: ",best_p4)

  pname=["mu_1","P_1","t01","ecos1","esin1",
        "mu_2","P_2","t02","ecos2","esin2",
        "mu_3","P_3","t03","ecos3","esin3",
        "mu_4","P_4","t04","ecos4","esin4"]
  mean_mp=[best_p4[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_p4[(iplanet-1)*5+4]^2 + best_p4[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  ecc_errs=[sqrt(err[(iplanet-1)*5+4]^2 + err[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]

  open(results,"w") do io
    println(io,"Global Fit Results. Per=[",p4in," - ",p4out,", length=",np4,"]")
    for i=1:nparam
      println(io,pname[i],": ",best_p4[i]," ± ",err[i])
    end
    println(io,"Retrieved Earth masses: ",mean_mp," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",mean_ecc," ± ",ecc_errs)
  end
  @save outfile p3 lprob_p3 best_p3 dp lprob_dp best_dp lprob_best_p3 p4 lprob_p4 best_p4 lprob_best_p4 ntrans nplanet tt0 tt ttmodel sigtt nphase
  return best_p4,lprob_best_p4 
end