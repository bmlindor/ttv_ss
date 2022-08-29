include("sim_times.jl")
include("misc.jl")
using TTVFaster,DelimitedFiles,JLD2,LsqFit,Statistics,DataFrames,CSV

function fit_planet5(filename::String,jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p3in::Float64,p3out::Float64,np3::Int,nphase::Int,p4in::Float64,p4out::Float64,np4::Int,p5in::Float64,p5out::Float64,np5::Int,obs::String)
  if obs=="fromEMB"
    fitfile = string("FITS/fromEMB/p5_fittest",sigma,"s",nyear,"yrs.jld2")
    results = string("results/fromEMB/p5_fittest",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/fromEMB/p5_gridtest",sigma,"s",nyear,"yrs.txt")
  elseif obs=="fromEV"
    fitfile = string("FITS/p5_fittest",sigma,"s",nyear,"yrs.jld2")
    results = string("results/p5_fittest",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/p5_gridtest",sigma,"s",nyear,"yrs.txt")
  end
   @assert isfile(filename)
  println(filename," loaded.")
  data1 = readdlm(filename,Float64,comments=true)
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
      param3 = [best_p2;param_tmp] #concatenate 2 planet model to 3 planet model params
      p3_cur = p3[j] #sets jupiter period to global value
      param1 = param3 .+ 100.0
      niter=0
      while maximum(abs.(param1 .- param3)) > tol && niter < 20
        param1 = param3
        fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:10];10^params[11];p3_cur;params[12:end]],jmax,true),tt0,tt,weight,param3)
        param3 = fit.param
        niter+=1
        # println("init_param: ",param3)
        # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,param3,tt,sigtt,true,p3_cur))
      end
      ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[param3[1:10];10^param3[11];p3_cur;param3[12:end]],jmax,true)
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
    # println("Period: ",p3[j]," log Prob: ",lprob_p3[j]," Param: ",vec(param_p3[1:nparam,j]))
  end
  println("Finished 3-planet fit w/ fixed period: ",p3best)

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,p3best)
  best_p3 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p3,jmax,true)
  lprob_best_p3= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  println("Finished global 3-planet fit in ",niter," iterations.")
  println("New 3-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p3,tt,sigtt,jmax,true))
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

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,p4best)
  best_p4 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p4,jmax,true)
  lprob_best_p4= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global 4-planet fit.")
  println("New 4-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p4,tt,sigtt,jmax,true))
  println("Maximum: ",lprob_best_p4," Param: ",best_p4)

 # Now,add a 5th planet:
  ntrans = [nt1,nt2,2,2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
  nplanet = 5
  nparam = 25
  # Grid of periods to search over:
  p5 = 10 .^ range(log10(p5in),stop=log10(p5out),length=np5)
  p5_cur =  1.88*365.25 ##29.44*365.25 
  lprob_p5 = zeros(np5)
  param_p5 = zeros(nparam,np5)
  lprob_best = -1e100 
  p5best = zeros(nparam)
  niter = 0
  for j=1:np5
    phase = p5[j]*range(0,stop=1,length=nphase) 
    lprob_phase = zeros(nphase) 
    lprob_p5[j] = -1e100
    for i=1:nphase
     # p5 param_names: mass ratio,phase,ecosw,esinw; uses same nphase as p3
      param_tmp = [log10(1e-7),phase[i],0.01,0.01]
      param5 = [best_p4[1:20];param_tmp]   
      p5_cur = p5[j]
      param1 = param5 .+ 100.0
      niter=0
      while maximum(abs.(param1 .- param5)) > tol && niter < 20
        param1 = param5
        # println("init_param: ",param5)
        # println(ntrans)
        fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:10];10^params[11];p5_cur;params[12:end]],jmax,true),tt0,tt,weight,param5)
        param5 = fit.param
        niter+=1 
      end
      ttmodel=ttv_wrapper(tt0,nplanet,ntrans,[param5[1:10];10^param5[11];p5_cur;param5[12:end]],jmax,true)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best
        lprob_best = lprob_phase[i]
        p5best = [fit.param[1:10];10^param5[11];p5_cur;fit.param[12:end]]
      end
      if lprob_phase[i] > lprob_p5[j] 
        lprob_p5[j] = lprob_phase[i]
        param_p5[1:nparam,j] =  [fit.param[1:10];10^param5[11];p5_cur;fit.param[12:end]]
      end
    end
    # println("Period: ",p5[j]," log Prob: ",lprob_p5[j]," Param: ",vec(param_p5[1:nparam,j]))
  end
  println("Finished 5-planet fit w/ fixed period: ",p5best," in ",niter," iterations")
  open(grid,"w") do io
    writedlm(io,zip(p5,lprob_p5))
  end
  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,p5best)
  cov=estimate_covar(fit)
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  best_p5 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p5,jmax,true)
  lprob_best_p5= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global 5-planet fit.")
  println("New 5-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p5,tt,sigtt,jmax,true))
  println("Maximum: ",lprob_best_p5," Param: ",best_p5)

  pname=["mu_1","P_1","t01","ecos1","esin1",
      "mu_2","P_2","t02","ecos2","esin2",
      "mu_3","P_3","t03","ecos3","esin3",
      "mu_4","P_4","t04","ecos4","esin4",
      "mu_5","P_5","t05","ecos5","esin5"]
  mean_mp=[best_p5[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_p5[(iplanet-1)*5+4]^2 + best_p5[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  ecc_errs=[sqrt(err[(iplanet-1)*5+4]^2 + err[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]

  open(results,"w") do io
    println(io,"Global Fit Results.",'\n',"per 5 range=[",p5in," - ",p5out,", length=",np5,"]")
    for i=1:nparam
      println(io,pname[i],": ",best_p5[i]," ± ",err[i])
    end
    println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",'\n',mean_ecc,'\n'," ± ",ecc_errs)
  end
  @save fitfile p4 lprob_p4 best_p4 lprob_best_p4 p5 lprob_p5 best_p5 lprob_best_p5 ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase p4in p4out np4 p5in p5out np5
  return best_p4, best_p5   
end

"""
    fit_planet5(jd1,sigma,nyear,tref,tol,p5in,p5out,np5,nphase,options)

# Arguments:
- `jd1::Float64`: starting Julian Ephemeris Date of observations.
- `sigma::Real`: fixed noised added to observations.
- `nyear::Real`: time span of observations.
- `tref::Real`: JED to subtract from transit times to aid fit of low mass planets.
- `tol::Real`: tolerance level of fit.
- `p5in::Float64`: starting period to perform seach for Saturn-like planet (in days)
- `p5out::Float64`: ending period to perform seach for Saturn-like planet (in days)
- `np5::Int`: number of periods to fit
- `nphase::Int`: number of phases to fit
- `options::Array{String}`:arg 1=source of observations for body 2 (EMB or EV); arg 2=whether grid is accurate or wide)
# Returns:
- `best_p3::Array{Float64}`: list of global best paramters for 4 planets given the observed transit times.
- `lprob_best_p3::Float64`: log probability of detecting 4 planets with the given properties.
"""
# If 4-planet fit already exists, can just do 5-planet search
function fit_planet5(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p5in::Float64,p5out::Float64,np5::Int,nphase::Int,options::Array{String},save_as_jld2::Bool=false)
	obs=options[1]; grid_type=options[2]
	if grid_type=="accurate"
	if obs=="fromEMB"
    infile = string("FITS/fromEMB/p4_fit",sigma,"s",nyear,"yrs.jld2")
    outfile = string("FITS/fromEMB/p5_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/fromEMB/p5_fit",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/fromEMB/p5_grid",sigma,"s",nyear,"yrs.csv")
  elseif obs=="fromEV"
    infile = string("FITS/p4_fit",sigma,"s",nyear,"yrs.jld2")
    outfile = string("FITS/p5_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/p5_fit",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/p5_grid",sigma,"s",nyear,"yrs.csv")
  end
	end
	if grid_type=="wide"
	if obs=="fromEMB"
    infile = string("FITS/fromEMB/p4_fit",sigma,"s",nyear,"yrs.jld2")
    outfile = string("FITS/fromEMB/widep5_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/fromEMB/widep5_fit",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/fromEMB/widep5_grid",sigma,"s",nyear,"yrs.csv")
  elseif obs=="fromEV"
    infile = string("FITS/p4_fit",sigma,"s",nyear,"yrs.jld2")
    outfile = string("FITS/widep5_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/widep5_fit",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/widep5_grid",sigma,"s",nyear,"yrs.csv")
  end
	end
  @assert isfile(infile)
  m = jldopen(String(infile),"r")
  tt0,tt,ttmodel,sigtt=m["tt0"],m["tt"],m["ttmodel"],m["sigtt"]
  nt1,nt2 = m["ntrans"][1],m["ntrans"][2]
  p4,lprob_p4=m["p4"],m["lprob_p4"]
  best_p4,lprob_best_p4=m["best_p4"],m["lprob_best_p4"]
  Nobs = sum([nt1,nt2])
  jmax=5
  jd2 = nyear*365.25 + jd1
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
  println(infile," loaded.")
  println("Previous model params: ",best_p4)

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
  niter = 0
  for j=1:np5
    phase = p5[j]*range(0,stop=1,length=nphase) 
    lprob_phase = zeros(nphase) 
    lprob_p5[j] = -1e100
    for i=1:nphase
     # p5 param_names: mass ratio,phase,ecosw,esinw; uses same nphase as p3
      param_tmp = [log10(1e-4),phase[i],0.01,0.01]
      param5 = [best_p4[1:20];param_tmp]   
      p5_cur = p5[j]
      param1 = param5 .+ 100.0
      niter=0
      while maximum(abs.(param1 .- param5)) > tol && niter < 20
        param1 = param5
        fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:20];10^params[21];p5_cur;params[22:end]],jmax,true),tt0,tt,weight,param5)
        param5 = fit.param 
        niter+=1
      end
      ttmodel=ttv_wrapper(tt0,nplanet,ntrans,[param5[1:20];10^param5[21];p5_cur;param5[22:end]],jmax,true)
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
  println("Finished 5-planet fit w/ fixed period: ",p5best," in ",niter," iterations")
	df=DataFrame(mu_1=param_p5[1,:],P_1=param_p5[2,:],t01=param_p5[3,:],ecos1=param_p5[4,:],esin1=param_p5[5,:],
								mu_2=param_p5[6,:],P_2=param_p5[7,:],t02=param_p5[8,:],ecos2=param_p5[9,:],esin2=param_p5[10,:],
								mu_3=param_p5[11,:],P_3=param_p5[12,:],t03=param_p5[13,:],ecos3=param_p5[14,:],esin3=param_p5[15,:],
								mu_4=param_p5[16,:],P_4=param_p5[17,:],t04=param_p5[18,:],ecos4=param_p5[19,:],esin4=param_p5[20,:],
								mu_5=param_p5[16,:],P_5=param_p5[17,:],t05=param_p5[18,:],ecos5=param_p5[19,:],esin5=param_p5[20,:],
								lprob=lprob_p5[:])
	CSV.write(grid,df)

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,p5best)
  cov=estimate_covar(fit) ;  best_p5 = fit.param 
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p5,jmax,true)
  lprob_best_p5= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global 5-planet fit.")
  println("New 5-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p5,tt,sigtt,jmax,true))
  println("Maximum: ",lprob_best_p5," Param: ",best_p5)

  pname=["mu_1","P_1","t01","ecos1","esin1",
      "mu_2","P_2","t02","ecos2","esin2",
      "mu_3","P_3","t03","ecos3","esin3",
      "mu_4","P_4","t04","ecos4","esin4",
      "mu_5","P_5","t05","ecos5","esin5"]
  mean_mp=[best_p5[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_p5[(iplanet-1)*5+4]^2 + best_p5[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  ecc_errs=[sqrt(err[(iplanet-1)*5+4]^2 + err[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]

  open(results,"w") do io
    println(io,"Global Fit Results.",'\n',"per 5 range=[",p5in," - ",p5out,", length=",np5,"]")
    for i=1:nparam
      println(io,pname[i],": ",best_p5[i]," ± ",err[i])
    end
    println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",'\n',mean_ecc,'\n'," ± ",ecc_errs)
  end
	if save_as_jld2
  @save outfile p4 lprob_p4 best_p4 lprob_best_p4 p5 lprob_p5 best_p5 lprob_best_p5 ntrans nplanet tt0 tt ttmodel sigtt nphase
	end
  return best_p5,lprob_best_p5   
end
