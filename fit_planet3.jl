# if !@isdefined(TTVFaster)
#     include("TTVFaster/src/TTVFaster.jl")
#     using Main.TTVFaster
# end
import Main.TTVFaster.ttv_wrapper
import Main.TTVFaster.chisquare
include("sim_times.jl")
include("CGS.jl")
using DelimitedFiles,JLD2,LsqFit,Statistics

function fit_planet3(filename::String,jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p3in::Float64,p3out::Float64,np3::Int,nphase::Int,obs::String)
  if obs=="fromEMB"
    fitfile = string("FITS/fromEMB/p3_fittest",sigma,"s",nyear,"yrs.jld2")
    results = string("results/fromEMB/p3_fittest",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/fromEMB/p3_gridtest",sigma,"s",nyear,"yrs.txt")
  elseif obs=="fromEV"
    fitfile = string("FITS/p3_fittest",sigma,"s",nyear,"yrs.jld2")
    results = string("results/p3_fittest",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/p3_gridtest",sigma,"s",nyear,"yrs.txt")
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
  ntrans = [nt1,nt2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
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
    # Loop over planet 3 phases:    
    for i=1:nphase
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
      if lprob_phase[i] > lprob_best 
      # check that best fit for period is better than global best fit
        lprob_best = lprob_phase[i]
        p3best = [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      end
      if lprob_phase[i] > lprob_p3[j] 
      # checks best fit over all phases of jupiter for this particular period
        lprob_p3[j] = lprob_phase[i]
        param_p3[1:nparam,j] =  [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      end
      # if j>1 && abs(lprob_p3[j] - lprob_p3[j-1])>5
      #   # Check that best fit for current period is close to that of previous period
      #   lprob_p3[j] = lprob_p3[j-1]
      #   param_p3[1:nparam,j] = [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      # end
    end
    # println("Period: ",p3[j]," log Prob: ",lprob_p3[j]," Param: ",vec(param_p3[1:nparam,j]))
  end
  println("Finished 3-planet fit w/ fixed period: ",p3best," in ",niter," iterations")
  # Make likelihood profile continuous???
      # for j=1:np3
      #   if abs(lprob_p3[j+1] - lprob_p3[j])
  # Rescale to set minimum chi-square equal to number of degrees of freedom
  #  = number of transits - number of model parameters (15):
  # plot(p3/365.25,exp.(-0.5*(chi_p3 .-minimum(chi_p3))))

  open(grid,"w") do io
    writedlm(io,zip(p3,lprob_p3))
  end
  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,p3best)
  cov=estimate_covar(fit)
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  best_p3 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p3,jmax,true)
  lprob_best_p3= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global 3-planet fit.")
  println("New 3-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p3,tt,sigtt,jmax,true))
  println("Maximum: ",lprob_best_p3," Param: ",best_p3)

  pname=["mu_1","P_1","t01","ecos1","esin1",
        "mu_2","P_2","t02","ecos2","esin2",
        "mu_3","P_3","t03","ecos3","esin3"]
  mean_mp=[best_p3[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_p3[(iplanet-1)*5+4]^2 + best_p3[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  ecc_errs=[sqrt(err[(iplanet-1)*5+4]^2 + err[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]

  open(results,"w") do io
    println(io,"Global Fit Results.",'\n',"per 3 range=[",p3in," - ",p3out,", length=",np3,"]")
    for i=1:nparam
      println(io,pname[i],": ",best_p3[i]," ± ",err[i])
    end
    println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",'\n',mean_ecc,'\n'," ± ",ecc_errs)
  end
  # plot_profile(p3,lprob_p3,p3_cur,50,[sigma,nyear],"firebrick","Jupiter")
  # title=string("IMAGES/likelihoods/",obs,"Jupiter",sigma,"secs",nyear,"yrs.png")
  # savefig(title)
  # Create files
  @save fitfile p3 lprob_p3 best_p3 lprob_best_p3 p3best ntrans nplanet tt0 tt ttmodel sigtt
  return lprob_best_p3,best_p3,lprob_p3,p3
end

"""
    fit_planet3(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,obs)

 If the 2-planet fit already exists, just do 3-planet fit.
# Arguments:
- `jd1::Float64`: starting Julian Ephemeris Date of observations.
- `sigma::Real`: fixed noised added to observations.
- `nyear::Real`: time span of observations.
- `tref::Real`: JED to subtract from transit times to aid fit of low mass planets.
- `tol::Real`: tolerance level of fit.
- `p3in::Float64`: starting period to perform seach for Jupiter-like planet (in days)
- `p3out::Float64`: ending period to perform seach for Jupiter-like planet (in days)
- `np3::Int`: number of periods to fit
- `nphase::Int`: number of phases to fit
- `obs::String`: source of observations for body 2 (EMB or EV).
# Returns:
- `best_p2::Vector{Float64}`: list of global best paramters for 3 planets given the observed transit times.
- `lprob_best_p2::Float64`: log probability of detecting 3 planets with the given properties.
"""
function fit_planet3(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p3in::Float64,p3out::Float64,np3::Int,nphase::Int,obs::String)
  if obs=="fromEMB"
    infile = string("FITS/fromEMB/p2_fit",sigma,"s",nyear,"yrs.jld2")
    outfile = string("FITS/fromEMB/p3_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/fromEMB/p3_fit",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/fromEMB/p3_grid",sigma,"s",nyear,"yrs.txt")
  elseif obs=="fromEV"
    infile = string("FITS/p2_fit",sigma,"s",nyear,"yrs.jld2")
    outfile = string("FITS/p3_fit",sigma,"s",nyear,"yrs.jld2")
    results = string("results/p3_fit",sigma,"s",nyear,"yrs.txt")
    grid = string("grid/p3_grid",sigma,"s",nyear,"yrs.txt")
  end
  @assert isfile(infile)
  p = jldopen(String(infile),"r")
  tt0,tt,ttmodel,sigtt=p["tt0"],p["tt"],p["ttmodel"],p["sigtt"]
  nt1,nt2 = p["ntrans"][1],p["ntrans"][2]
  best_p2=p["best_p2"]
  Nobs = sum([nt1,nt2])
  jmax=5
  jd2 = nyear*365.25 + jd1
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
  println(infile," loaded.")
  println("Previous model params: ",best_p2)
  
  # Now,let's add the 3rd planet:
  ntrans = [nt1,nt2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
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
    # Loop over planet 3 phases:
    for i=1:nphase 
     # p3 param_names: mass ratio,phase,ecosw,esinw
      param_tmp = [log10(1e-3),phase[i],0.01,0.01] 
      param3 = [best_p2;param_tmp] 
      p3_cur = p3[j] 
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
        p3best = [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      end
      if lprob_phase[i] > lprob_p3[j] # checks best fit over all phases of jupiter for this particular period
        lprob_p3[j] = lprob_phase[i]
        param_p3[1:nparam,j] = [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      end
    end
    # println("Period: ",p3[j]," log Prob: ",lprob_p3[j]," Param: ",vec(param_p3[1:nparam,j]))
  end
  println("Finished 3-planet fit w/ fixed period: ",p3best," in ",niter," iterations")
	open(grid,"w") do io
	  writedlm(io,zip(p3,lprob_p3))
	end
  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,p3best)
  cov=estimate_covar(fit)
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  best_p3 = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p3,jmax,true)
  lprob_best_p3= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global 3-planet fit.")
  println("New 3-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_p3,tt,sigtt,jmax,true))
  println("Maximum: ",lprob_best_p3," Param: ",best_p3)

  pname=["mu_1","P_1","t01","ecos1","esin1",
        "mu_2","P_2","t02","ecos2","esin2",
        "mu_3","P_3","t03","ecos3","esin3"]
  mean_mp=[best_p3[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_p3[(iplanet-1)*5+4]^2 + best_p3[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  ecc_errs=[sqrt(err[(iplanet-1)*5+4]^2 + err[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]

  open(results,"w") do io
    println(io,"Global Fit Results.",'\n',"per 3 range=[",p3in," - ",p3out,", length=",np3,"]")
    for i=1:nparam
      println(io,pname[i],": ",best_p3[i]," ± ",err[i])
    end
    println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",'\n',mean_ecc,'\n'," ± ",ecc_errs)
  end
  @save outfile p3 lprob_p3 best_p3 lprob_best_p3 ntrans nplanet tt0 tt ttmodel sigtt nphase
  return best_p3,lprob_best_p3
end
