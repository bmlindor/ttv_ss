include("sim_times.jl")
include("misc.jl")
include("CGS.jl")
using TTVFaster,DelimitedFiles,JLD2,LsqFit,Statistics,DataFrames,CSV

"""
    fit_moon(jd1,sigma,nyear,tref,tol,p3in,p3out,np3,nphase,options)
# Arguments:
- `jd1::Float64`: starting Julian Ephemeris Date of observations.
- `sigma::Real`: fixed noised added to observations.
- `nyear::Real`: time span of observations.
- `tref::Real`: JED to subtract from transit times to aid fit of low mass planets.
- `tol::Real`: tolerance level of fit.
- `dpin::Float64`: starting deltaphi phase to perform seach for Moon
- `dpout::Float64`: ending deltaphi phase to perform seach for Moon
- `ndp::Int`: number of phases to fit
- `options::Array{String}`:arg 1=whether grid is accurate or wide
- `nplanets::Real`: number of planets to consider
# Returns:
- `best_dp::Array{Float64}`: list of global best paramters for nplanets+moon given the observed transit times.
- `lprob_best_dp::Float64`: log probability of detecting nplanets+moon with the given properties.
"""
# If planet fit already exists, can just do moon search
function fit_moon(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,dpin::Float64,dpout::Float64,ndp::Int,options::Array{String},save_as_jld2::Bool=false,nplanets::Real=3)
	grid_type_nplanet=options[1]	
  infile = string("../FITS/p",nplanets,"_fit",sigma,"s",nyear,"yrs.jld2")
  outfile = string("../FITS/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
  results = string("../results/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.txt")
  grid = string("../grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
  @assert isfile(infile)
  m = jldopen(String(infile),"r")
  tt0,tt,ttmodel,sigtt=m["tt0"],m["tt"],m["ttmodel"],m["sigtt"]
  nt1,nt2 = m["ntrans"][1],m["ntrans"][2]
  nplanet,ntrans = m["nplanet"],m["ntrans"]
  best_par = m[string("best_p",nplanets)]
  lprob_best_par = m[string("lprob_best_p",nplanets)]
  if nplanets>2
    per, lprob_per = m[string("p",nplanets)],m[string("lprob_p",nplanets)]
  end
  Nobs = sum([nt1,nt2])
  jmax=5
  jd2 = nyear*365.25 + jd1
  weight = ones(nt1+nt2)./ sigtt.^2
  println(infile," loaded.")
  println("Previous model params: ",best_par)

  # Now,search for Moon:
  nparam = length(best_par)+3
  dp_cur = 2.312
  dp = range(dpin,stop=dpout,length=ndp)
  lprob_dp = zeros(ndp)
  param_dp = zeros(nparam,ndp)
  lprob_best = -1e100 
  dpbest = zeros(nparam)
  niter = 0
  for j=1:ndp
    lprob_dp[j] = -1e100 
    # lunar params: tmax sin(phi_0) ,tmax cos(phi_0) ,deltaphi 
    param_tmp = [0.01,0.01,dp[j]] 
    param5 = [best_par;param_tmp]
    dp_cur = dp[j]
    param1 = param5 .+ 100.0
		niter=0
    while maximum(abs.(param1 .- param5)) > tol && niter < 20
      param1 = param5
      fit = curve_fit((tt0,param5) -> ttv_wrapper(tt0,nplanet,ntrans,param5,jmax,false),tt0,tt,weight,param5)
      param5 = fit.param
      niter+=1 
    end
    ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[fit.param[1:end-1];dp_cur],jmax,false)
    lprob_dp[j]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
    if lprob_dp[j] > lprob_best 
      lprob_best = lprob_dp[j]
      dpbest = [fit.param[1:end-1];dp_cur]
    end
    param_dp[1:nparam,j] = [fit.param[1:end-1];dp_cur]
    # println("deltaphi: ",deltaphi[j]," log Prob: ",lprob_dp[j]," Param: ",vec(param_dp[1:nparam,j]))
  end
  println("Finished lunar search: ",dpbest," in ",niter," iterations")
	df=DataFrame(mu_1=param_dp[1,:],P_1=param_dp[2,:],t01=param_dp[3,:],ecos1=param_dp[4,:],esin1=param_dp[5,:],
								mu_2=param_dp[6,:],P_2=param_dp[7,:],t02=param_dp[8,:],ecos2=param_dp[9,:],esin2=param_dp[10,:],
								mu_3=param_dp[11,:],P_3=param_dp[12,:],t03=param_dp[13,:],ecos3=param_dp[14,:],esin3=param_dp[15,:],
								tcos=param_dp[end-2,:],tsin=param_dp[end-1,:],dphi=param_dp[end,:],
								lprob=lprob_dp[:])
	CSV.write(grid,df) 

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,false),tt0,tt,weight,dpbest)
  cov=estimate_covar(fit) ;  best_dp = fit.param 
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_dp,jmax,false)
  lprob_best_dp = (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global moon fit.")
  chi2=chisquare(tt0,nplanet,ntrans,best_dp,tt,sigtt,jmax,false)
  println("Lunar chi-square: ",chi2)
  println("Maximum: ",lprob_best_dp," Param: ",best_dp)

  pname=["mu_1","P_1","t01","ecos1","esin1",
        "mu_2","P_2","t02","ecos2","esin2",
        "mu_3","P_3","t03","ecos3","esin3",
				"tcos","tsin","dphi"]
  mean_mp=[best_dp[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_dp[(iplanet-1)*5+4]^2 + best_dp[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  
  open(results,"w") do io
    println(io,"Global Fit Results.",'\n',"chi^2: ",chi2,'\n',"Δϕ range=[",dpin," - ",dpout,", length=",ndp,"]")
    for i=1:nparam
      println(io,pname[i],": ",best_dp[i]," ± ",err[i])
    end
    println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",'\n',mean_ecc)
  end
	if save_as_jld2
  @save outfile dp lprob_dp best_dp lprob_best_dp ntrans nplanet tt0 tt ttmodel sigtt
	end
  return best_dp,lprob_best_dp
end

"""
    fit_moon(jd1,sigma,nyear,tref,tol,p4in,p4out,np4,nphase,options)

# Arguments:
- `jd1::Float64`: starting Julian Ephemeris Date of observations.
- `sigma::Real`: fixed noised added to observations.
- `nyear::Real`: time span of observations.
- `tref::Real`: JED to subtract from transit times to aid fit of low mass planets.
- `tol::Real`: tolerance level of fit.
- `p4in::Float64`: starting period to perform seach for Mars-like planet (in days)
- `p4out::Float64`: ending period to perform seach for Mars-like planet (in days)
- `np4::Int`: number of periods to fit
- `nphase::Int`: number of phases to fit
- `options::Array{String}`:arg 1=whether grid is accurate or wide
# Returns:
- `best_p4::Array{Float64}`: list of global best paramters for 4 planets given the observed transit times.
- `lprob_best_p4::Float64`: log probability of detecting 4 planets with the given properties.
"""
# If 3-planet fit with moon already exists, can do 4-planet search
function fit_moon(jd1::Float64,sigma::Real,nyear::Real,tref::Real,tol::Real,p4in::Float64,p4out::Float64,np4::Int,nphase::Int,options::Array{String},save_as_jld2::Bool=false)
	grid_type_nplanet=options[1]
  infile = string("../FITS/p3moon_fit",sigma,"s",nyear,"yrs.jld2")
  outfile = string("../FITS/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
  results = string("../results/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.txt")
  grid = string("../grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
  @assert isfile(infile)
  m = jldopen(String(infile),"r")
  tt0,tt,ttmodel,sigtt=m["tt0"],m["tt"],m["ttmodel"],m["sigtt"]
  nt1,nt2 = m["ntrans"][1],m["ntrans"][2]
  dp = m["dp"] ;  lprob_dp=m["lprob_dp"]
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
      param_tmp = [log10(1e-7),phase[i],0.01,0.01]
      # Mars' period is shorter than Jupiter's, so need to keep sorted for now
      param4 = [best_dp[1:10];param_tmp;best_dp[11:end]]   
      p4_cur = p4[j]
      param1 = param4 .+ 100.0
      niter=0
      while maximum(abs.(param1 .- param4)) > tol && niter <20
        param1 = param4
        fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:10];10^params[11];p4_cur;params[12:end]],jmax,false),tt0,tt,weight,param4)
        param4 = fit.param 
        niter+=1
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
	df=DataFrame(mu_1=param_p4[1,:],P_1=param_p4[2,:],t01=param_p4[3,:],ecos1=param_p4[4,:],esin1=param_p4[5,:],
								mu_2=param_p4[6,:],P_2=param_p4[7,:],t02=param_p4[8,:],ecos2=param_p4[9,:],esin2=param_p4[10,:],
								mu_3=param_p4[11,:],P_3=param_p4[12,:],t03=param_p4[13,:],ecos3=param_p4[14,:],esin3=param_p4[15,:],
								mu_4=param_p4[16,:],P_4=param_p4[17,:],t04=param_p4[18,:],ecos4=param_p4[19,:],esin4=param_p4[20,:],
								tcos=param_p4[end-2,:],tsin=param_p4[end-1,:],dphi=param_p4[end,:],
								lprob=lprob_p4[:])
	CSV.write(grid,df)

  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,false),tt0,tt,weight,p4best)
  cov=estimate_covar(fit) ;  best_p4 = fit.param 
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_p4,jmax,false)
  lprob_best_p4= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global 4-planet fit.")
  chi2=chisquare(tt0,nplanet,ntrans,best_p4,tt,sigtt,jmax,false)
	println("New 4-planet chi-square: ",chi2)
  println("Maximum: ",lprob_best_p4," Param: ",best_p4)

  pname=["mu_1","P_1","t01","ecos1","esin1",
        "mu_2","P_2","t02","ecos2","esin2",
        "mu_3","P_3","t03","ecos3","esin3",
        "mu_4","P_4","t04","ecos4","esin4",
				"tsin","tcos","dphi"]
  mean_mp=[best_p4[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_p4[(iplanet-1)*5+4]^2 + best_p4[(iplanet-1)*5+5]^2) for iplanet=1:nplanet]

  open(results,"w") do io
    println(io,"Global Fit Results.",'\n',"chi^2: ",chi2,'\n',"per 4 range=[",p4in," - ",p4out,", length=",np4,"]")
    for i=1:nparam
      println(io,pname[i],": ",best_p4[i]," ± ",err[i])
    end
    println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",'\n',mean_ecc)
  end
	if save_as_jld2
  @save outfile dp lprob_dp best_dp lprob_best_dp p4 lprob_p4 best_p4 lprob_best_p4 ntrans nplanet tt0 tt ttmodel sigtt nphase
	end
  return best_p4,lprob_best_p4 
end
