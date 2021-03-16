if !@isdefined(TTVFaster)
    include("TTVFaster/src/TTVFaster.jl")
    using Main.TTVFaster
end
import Main.TTVFaster.ttv_wrapper
import Main.TTVFaster.chisquare
include("regress.jl")
using DelimitedFiles,JLD2,Optim,LsqFit,Statistics

function fit_planet3(filename::String,
  jd1::Float64,nyear::Float64,jdsize::Int64,
  p3in::Float64,p3out::Float64,np3::Int,nphase::Int,
  addnoise::Bool=false,sigma::Float64=0.0,EMB::Bool=true)
  
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
  p3_cur = 11.86*365.25 #jupiter period in days,initial value
  #res = optimize(chisquare2,param,method = :l_bfgs,iterations = 21)
  ntrans = [nt1,nt2]
  Nobs = sum(ntrans)
  nplanet = 2
  # create initial simplex? need function for this?
  # result = optimize(f0,xcurr,NelderMead(initial_simplex=MySimplexer(),show_trace=true,iterations=1))
  println("Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax))
  # res = optimize(params -> chisquare(nplanet,ntrans,params,tt,sigtt),init_param) 
  # init_param = res.minimizer
  param1 = init_param .+ 100.0
  while maximum(abs.(param1 .- init_param)) > 1e-5
    param1 = init_param
    res = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax),tt0,tt,weight,init_param)
    init_param = res.param
    # println("init_param: ",init_param)
    # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt))
  end
  # res = optimize(params -> chisquare(tt0,nplanet,ntrans,params,tt,sigtt),init_param) 
  # init_param = res.minimizer
  # fit2 = curve_fit(ttv_wrapper2,tt0,tt,weight,param; show_trace=true)
  println("Finished 2-planet fit: ",init_param)

  # Now,let's add the 3rd planet:
  ntrans = [nt1,nt2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
  nplanet = 3
  #p3 = 11.86*365.25
  # Grid of periods to search over:
  p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
  lprob_p3 = zeros(np3)
  nparam = 15
  param_p3 = zeros(nparam,np3)
  lprob_best = -1e100 #global best fit
  pbest = zeros(nparam)
  # Shifting to simulated observation range to search over period grid
  offset = (jd1 + jd2)/2 
  for j=1:np3
    phase = p3[j]*range(0,stop=1,length=nphase) .+ offset 
    lprob_phase = zeros(nphase)
    lprob_p3[j] = -1e100
    for i=1:nphase #loops over jupiter phases
      param_tmp = [1e-3,phase[i],0.01,0.01] # jupiter params: mass ratio,phase,ecosw,esinw
      param3 = [init_param;param_tmp] #concatenate 2 planet model to 3 planet model params
      p3_cur = p3[j] #sets jupiter period to global value
      # fit = curve_fit(ttv_wrapper_fixp3,tt0,tt,weight,param3) #optimizes fit w/ 3 planet model
      # fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,true,p3_cur),tt0,tt,weight,param3) 
      # param3 = fit.param
      param1 = param3 .+ 100.0
      while maximum(abs.(param1 .- param3)) > 1e-5
        param1 = param3
        fit = curve_fit((tt0,param3) -> ttv_wrapper(tt0,nplanet,ntrans,[param3[1:11];p3_cur;param3[12:end]],jmax,true),tt0,tt,weight,param3)
        param3 = fit.param
        # println("init_param: ",param3)
        # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,param3,tt,sigtt,true,p3_cur))
      end
      ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[param3[1:11];p3_cur;param3[12:end]],jmax,true)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best # check that best fit for period is better than global best fit
        lprob_best = lprob_phase[i]
        pbest = [fit.param[1:11];p3_cur;fit.param[12:14]]
      end
      if lprob_phase[i] > lprob_p3[j] # checks best fit over all phases of jupiter for this particular period
        lprob_p3[j] = lprob_phase[i]
        param_p3[1:nparam,j] =  [fit.param[1:11];p3_cur;fit.param[12:14]]
      end
    end
    println("Period: ",p3[j]," chi: ",lprob_p3[j]," Param: ",vec(param_p3[1:nparam,j]))
  end
  println("Finished 3-planet fit w/ fixed period: ",pbest)

  #ttmodel=ttv_wrapper3(tt0,param3)
  #res = optimize(chisquare3,param3,method = :l_bfgs,iterations = 21)
  #  res = optimize(chisquare3,param3,method = :l_bfgs)
  #  ttmodel=ttv_wrapper3(tt0,param3)

  # fit = curve_fit(ttv_wrapper3,tt0,tt,weight,pbest)
  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax),tt0,tt,weight,pbest)
  # ttmodel=ttv_wrapper3(tt0,pbest)
  pbest_global = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,pbest_global,jmax)
  lprob_best= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  sigsys2 = 1e-6

  println("Finished global 3-planet fit.")
  println("Maximum: ",lprob_best," Param: ",pbest_global)


  pname = ["mu_1","P_1","t01","e1 cos(om1)","e1 sin(om1)",
        "mu_2","P_2","t02","e2 cos(om2)","e2 sin(om2)",
        "mu_3","P_3","t03","e3 cos(om3)","e3 sin(om3)"]

  results = string("OUTPUTS/p3_fit",sigma,"s",nyear,"yrs.txt")
  open(results,"w") do io
    for i=1:nparam
      println(io,pname[i],": ",pbest_global[i])
    end
  end
  fitfile = string("FITS/p3_fit",sigma,"s",nyear,"yrs.jld2")
  @save fitfile param_p3 lprob_p3 lprob_best pbest_global ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
  # writedlm(results,pbest_global)
    return lprob_best,pbest_global
end
