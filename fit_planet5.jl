# if !@isdefined(TTVFaster)
#     include("TTVFaster/src/TTVFaster.jl")
#     using Main.TTVFaster
# end
# import Main.TTVFaster.ttv_wrapper
# import Main.TTVFaster.chisquare
# include("regress.jl")
# using DelimitedFiles,JLD2,Optim,LsqFit,Statistics

function fit_planet5(filename::String,
  jd1::Float64,sigma::Float64,nyear::Float64,
  p3in::Float64,p3out::Float64,np3::Int,nphase::Int,
  p4in::Float64,p4out::Float64,np4::Int,
  p5in::Float64,p5out::Float64,np5::Int,
  addnoise::Bool=true,EM::Bool=true)
  jd2 = nyear*365.25 + jd1
  data1 = readdlm(filename)
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
  while maximum(abs.(param1 .- init_param)) > 1e-5
    param1 = init_param
    res = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,init_param)
    init_param = res.param
    # println("init_param: ",init_param)
    # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt))
  end
  # res = optimize(params -> chisquare(tt0,nplanet,ntrans,params,tt,sigtt),init_param) 
  # init_param = res.minimizer
  # fit2 = curve_fit(ttv_wrapper2,tt0,tt,weight,param; show_trace=true)
  println("Finished 2-planet fit: ",init_param)
  println("New p2 chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,EM))

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
  # Shifting to simulated observation range to search over period grid
  offset = (jd1 + jd2)/2 
  for j=1:np3
    phase = p3[j]*range(0,stop=1,length=nphase) .+ offset 
    lprob_phase = zeros(nphase)
    lprob_p3[j] = -1e100
    for i=1:nphase #loops over jupiter phases
     # p3 param_names: mass ratio,phase,ecosw,esinw
      param_tmp = [log10(1e-3),phase[i],0.01,0.01] 
      param3 = [init_param;param_tmp] #concatenate 2 planet model to 3 planet model params
      p3_cur = p3[j] #sets jupiter period to global value
      param1 = param3 .+ 100.0
      while maximum(abs.(param1 .- param3)) > 1e-5
        param1 = param3
        fit = curve_fit((tt0,param3) -> ttv_wrapper(tt0,nplanet,ntrans,[param3[1:10];10^param3[11];p3_cur;param3[12:end]],jmax,EM),tt0,tt,weight,param3)
        param3 = fit.param
        # println("init_param: ",param3)
        # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,param3,tt,sigtt,true,p3_cur))
      end
      ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[param3[1:10];10^param3[11];p3_cur;param3[12:end]],jmax,EM)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best # check that best fit for period is better than global best fit
        lprob_best = lprob_phase[i]
        p3best = [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:14]]
      end
      if lprob_phase[i] > lprob_p3[j] # checks best fit over all phases of jupiter for this particular period
        lprob_p3[j] = lprob_phase[i]
        param_p3[1:nparam,j] =  [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:14]]
      end
    end
    println("Period: ",p3[j]," log Prob: ",lprob_p3[j]," Param: ",vec(param_p3[1:nparam,j]))
  end
  println("Finished 3-planet fit w/ fixed period: ",p3best)

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
  offset = (jd1 + jd2)/2 
  for j=1:np4
    phase = p4[j]*range(0,stop=1,length=nphase) .+ offset 
    lprob_phase = zeros(nphase) 
    lprob_p4[j] = -1e100
    for i=1:nphase
     # p4 param_names: mass ratio,phase,ecosw,esinw; uses same nphase as p3
      param_tmp = [log10(1e-7),phase[i],0.01,0.01]
      # Mars' period is shorter than Jupiter's, so need to keep sorted for now
      param4 = [best_p3[1:10];param_tmp;best_p3[11:15]]   
      p4_cur = p4[j]
      param1 = param4 .+ 100.0
      while maximum(abs.(param1 .- param4)) > 1e-5
        param1 = param4
        fit = curve_fit((tt0,param4) -> ttv_wrapper(tt0,nplanet,ntrans,[param4[1:10];10^param4[11];p4_cur;param4[12:end]],jmax,EM),tt0,tt,weight,param4)
        param4 = fit.param 
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
    println("Period: ",p4[j]," log Prob: ",lprob_p4[j]," Param: ",vec(param_p4[1:nparam,j]))
  end
  println("Finished 4-planet fit w/ fixed period: ",p4best)

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,p4best)
  best_p4 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p4,jmax,EM)
  lprob_best_p4= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  println("Finished global 4-planet fit.")
  println("New 4-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p4,tt,sigtt,jmax,EM))
  println("Maximum: ",lprob_best_p4," Param: ",best_p4)

 # Now,add a 5th planet:
  ntrans = [nt1,nt2,2,2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
  nplanet = 5
  nparam = 25
  # Grid of periods to search over:
  p5 = 10 .^ range(log10(p5in),stop=log10(p5out),length=np5)
  p5_cur =  29.44*365.25 
  lprob_p5 = zeros(np5)
  param_p5 = zeros(nparam,np5)
  lprob_best = -1e100 
  p5best = zeros(nparam)
  offset = (jd1 + jd2)/2 
  for j=1:np5
    phase = p5[j]*range(0,stop=1,length=nphase) .+ offset 
    lprob_phase = zeros(nphase) 
    lprob_p5[j] = -1e100
    for i=1:nphase
     # p5 param_names: mass ratio,phase,ecosw,esinw; uses same nphase as p3
      param_tmp = [log10(1e-4),phase[i],0.01,0.01]
      param5 = [best_p4[1:20];param_tmp]   
      p5_cur = p5[j]
      param1 = param5 .+ 100.0
      while maximum(abs.(param1 .- param5)) > 1e-5
        param1 = param5
        # println("init_param: ",param5)
        # println(ntrans)
        fit = curve_fit((tt0,param5) -> ttv_wrapper(tt0,nplanet,ntrans,[param5[1:20];10^param5[21];p5_cur;param5[22:end]],jmax,EM),tt0,tt,weight,param5)
        param5 = fit.param 
      end
      ttmodel=ttv_wrapper(tt0,nplanet,ntrans,[param5[1:20];10^param5[21];p5_cur;param5[22:end]],jmax,EM)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best
        lprob_best = lprob_phase[i]
        p5best = [fit.param[1:20];10^param5[21];p5_cur;fit.param[22:end]]
      end
      if lprob_phase[i] > lprob_p5[j] 
        lprob_p5[j] = lprob_phase[i]
        param_p5[1:nparam,j] =  [fit.param[1:20];10^param5[21];p5_cur;fit.param[22:end]]
      end
    end
    println("Period: ",p5[j]," log Prob: ",lprob_p5[j]," Param: ",vec(param_p5[1:nparam,j]))
  end
  println("Finished 5-planet fit w/ fixed period: ",p5best)

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,p5best)
  best_p5 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p5,jmax,EM)
  lprob_best_p5= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  println("Finished global 5-planet fit.")
  println("New 5-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p5,tt,sigtt,jmax,EM))
  println("Maximum: ",lprob_best_p5," Param: ",best_p5)
  if EM
    fitfile = string("FITS/fromEMB/p5_fit",sigma,"s",nyear,"yrs.jld2")
  else
    fitfile = string("FITS/p5_fit",sigma,"s",nyear,"yrs.jld2")
  end  
  @save fitfile p4 lprob_p4 best_p4 lprob_best_p4 p5 lprob_p5 best_p5 lprob_best_p5 ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase p4in p4out np4 p5in p5out np5
  return best_p4, best_p5   
end
# If 4-planet fit already exists, can just do 5-planet search but what about after moon?
function fit_planet5(jd1::Float64,sigma::Float64,nyear::Float64,
  p5in::Float64,p5out::Float64,np5::Int,nphase::Int
  ,EM::Bool=true)
  if EM
    infile = string("FITS/fromEMB/p4_fit",sigma,"s",nyear,"yrs.jld2")
  else
    infile = string("FITS/p4_fit",sigma,"s",nyear,"yrs.jld2")
  end
  m = jldopen(String(infile),"r")
  tt0,tt,ttmodel,sigtt=m["tt0"],m["tt"],m["ttmodel"],m["sigtt"]
  nt1,nt2 = m["ntrans"][1],m["ntrans"][2]
  p4,lprob_p4=m["p4"],m["lprob_p4"]
  best_p4,lprob_best_p4=m["best_p4"],m["lprob_best_p4"]
  Nobs = sum([nt1,nt2])
  jmax=5
  jd2 = nyear*365.25 + jd1
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
 # Now,add a 5th planet:
  ntrans = [nt1,nt2,2,2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
  nplanet = 5
  nparam = 25
  println("Planet 4 fit loaded.")
  # Grid of periods to search over:
  p5 = 10 .^ range(log10(p5in),stop=log10(p5out),length=np5)
  p5_cur =  29.44*365.25 
  lprob_p5 = zeros(np5)
  param_p5 = zeros(nparam,np5)
  lprob_best = -1e100 
  p5best = zeros(nparam)
  offset = (jd1 + jd2)/2 
  for j=1:np5
    phase = p5[j]*range(0,stop=1,length=nphase) .+ offset 
    lprob_phase = zeros(nphase) 
    lprob_p5[j] = -1e100
    for i=1:nphase
     # p5 param_names: mass ratio,phase,ecosw,esinw; uses same nphase as p3
      param_tmp = [log10(1e-4),phase[i],0.01,0.01]
      param5 = [best_p4[1:20];param_tmp]   
      p5_cur = p5[j]
      param1 = param5 .+ 100.0
      while maximum(abs.(param1 .- param5)) > 1e-5
        param1 = param5
        # println("init_param: ",param5)
        # println(ntrans)
        fit = curve_fit((tt0,param5) -> ttv_wrapper(tt0,nplanet,ntrans,[param5[1:20];10^param5[21];p5_cur;param5[22:end]],jmax,EM),tt0,tt,weight,param5)
        param5 = fit.param 
      end
      ttmodel=ttv_wrapper(tt0,nplanet,ntrans,[param5[1:20];10^param5[21];p5_cur;param5[22:end]],jmax,EM)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best
        lprob_best = lprob_phase[i]
        p5best = [fit.param[1:20];10^param5[21];p5_cur;fit.param[22:end]]
      end
      if lprob_phase[i] > lprob_p5[j] 
        lprob_p5[j] = lprob_phase[i]
        param_p5[1:nparam,j] =  [fit.param[1:20];10^param5[21];p5_cur;fit.param[22:end]]
      end
    end
    # println("Period: ",p5[j]," log Prob: ",lprob_p5[j]," Param: ",vec(param_p5[1:nparam,j]))
  end
  println("Finished 5-planet fit w/ fixed period: ",p5best)

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,p5best)
  best_p5 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p5,jmax,EM)
  lprob_best_p5= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  println("Finished global 5-planet fit.")
  println("New 5-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p5,tt,sigtt,jmax,EM))
  println("Maximum: ",lprob_best_p5," Param: ",best_p5)
  if EM
    fitfile = string("FITS/fromEMB/p5_fit",sigma,"s",nyear,"yrs.jld2")
  else
    fitfile = string("FITS/p5_fit",sigma,"s",nyear,"yrs.jld2")
  end  
  @save fitfile p4 lprob_p4 best_p4 lprob_best_p4 p5 lprob_p5 best_p5 lprob_best_p5 ntrans nplanet tt0 tt ttmodel sigtt 
  return best_p4, best_p5   
end