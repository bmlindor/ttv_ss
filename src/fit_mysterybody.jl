# Julia v1.1
if !@isdefined(TTVFaster)
    include("TTVFaster/src/TTVFaster.jl")   
end
# using Main.TTVFaster
import Main.TTVFaster.ttv_wrapper
import Main.TTVFaster.chisquare
include("regress.jl")
include("misc.jl")
include("CGS.jl")
include("plot_likelihood.jl")
using DelimitedFiles,JLD2,LsqFit,Statistics,Profile,PyPlot
# using PyPlot,Unitful,UnitfulAstro,LinearAlgebra

jd1=2.4332825e6
tref = 2430000
tol = 1e-5
mratio = 1e-7
per_guess = 1.88*365.25
per_in, per_out, nper, nphase = 1.6*365.25, 3*365.25, 10, 36
jmax = 5
planet="Mars"
datafile="INPUTS/tt_30.0sEMB30.0yrs.txt"

function fit_mysteryplanet() #fit_mysteryplanet(datafile::String,jd1::Float64,tref::Real,tol::Real,obs::String)#::Int,mratio::Float64,per_guess::Float64,per_in::Float64,per_out::Float64,nper::Int,nphase::Int,
  #outfile = string("FITS/",obs,"/mystery_",planet,"_fit",jd1,"JED.jld2")
  outfile = string("FITS/mystery_",planet,"_fit",jd1,".txt")
  @assert isfile(datafile)
  println(datafile," loaded.")
  data1 = readdlm(datafile,Float64)
  tt = data1[:,3] .- tref
  sigtt = data1[:,4]
  nplanet_cond = 2
	nparam=10
  ntrans = zeros(nplanet_cond)
  ntrans=[sum(data1[:,1] .== i) for i=1:nplanet_cond]
  Nobs = sum(ntrans)
  # Actual transit times:
  nt1,nt2 = ntrans
  tt1,sigtt1 = tt[1:nt1], sigtt[1:nt1]
  tt2,sigtt2 = tt[nt1+1:nt1+nt2], sigtt[nt1+1:nt1+nt2]
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
  # Estimate the periods of each conditioned planet:
  p1est = median(tt1[2:end] - tt1[1:end-1])
  p2est = median(tt2[2:end] - tt2[1:end-1])
  # Okay,let's do a linear fit to the transit times:
  coeff1,covcoeff1 = find_coeffs(tt1,p1est,sigtt1)
  coeff2,covcoeff2 = find_coeffs(tt2,p2est,sigtt2)
  # Best-fit linear transit times without TTVs:
  t01 = coeff1[1]; per1 = coeff1[2]
  t02 = coeff2[1]; per2 = coeff2[2]
  t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) 
  t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
  tt0 = [t1;t2]

  # Okay,now let's do a 2-planet fit:
    # param_names = mass ratio,period,initial transit time,e*cos(omega),e*sin(omega)
	init_param_guess = [3e-6,per1,t01,0.01,0.01,
                  3e-6,per2,t02,0.01,0.01] 
    # println("Initial parameters: ",init_param)
	function global_fit(tt0,nplanet,ntrans,init_param,jmax,EM) 
	  # Perform fit with best params, and calculate covariances for parameters: 
  fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,init_param)
  covar=estimate_covar(fit)
  best_global = fit.param ;   nparam=length(best_global)
  err=[sqrt(covar[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_global,jmax,EM)
  lprob_best_global= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global fit.")
  println("New chi-square: ",chisquare(tt0,nplanet,ntrans,best_global,tt,sigtt,jmax,EM))
  println("Maximum: ",lprob_best_global,'\n'," Param: ",best_global)
  return best_global,err
  end

  function fit_cond_planets(tt0,tt,sigtt,jmax,nplanet,init_param)
    # Set up data structure to hold planet properties,passed to TTVFaster
    data=init_param
    p1=TTVFaster.Planet_plane_hk(data[1],data[2],data[3],data[4],data[ 5])
    p2=TTVFaster.Planet_plane_hk(data[6],data[7],data[8],data[9],data[10])
    # assuming no transits are skipped/duplicated ########### To Change
    time1 = collect(p1.trans0 .+ range(0,stop=nt1-1,length=nt1) .* p1.period)
    time2 = collect(p2.trans0 .+ range(0,stop=nt2-1,length=nt2) .* p2.period)
    # Initialize the computation of the Laplace coefficients:
    ttv1 = zeros(nt1)
    ttv2 = zeros(nt2)
    # Need first call to TTVFaster,without optimizing
    dummy=TTVFaster.compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2)

    # Now,optimize 2-planet fit
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

		best_p2,err=global_fit(tt0,nplanet,ntrans,init_param_guess,jmax,true)
    return best_p2,err
  end
  @time best_p2,err = fit_cond_planets(tt0,tt,sigtt,jmax,nplanet_cond,init_param_guess)
  # Now,let's add the 3rd planet:
  ntrans = [ntrans;2] #requires at least 2 transits for each planet (even if it doesnt transit)
  nplanet = nplanet_cond + 1
  nparam = length(best_p2) + 5
  # Grid of periods to search over:
  per = 10 .^ range(log10(per_in),stop=log10(per_out),length=nper)
  lprob_per = zeros(nper)
  per_cur = per_guess #jupiter period in days,initial value
  param_per = zeros(nparam,nper)
  lprob_best = -1e100 #global best fit
  perbest = zeros(nparam)
  niter = 0
  for j=1:nper
    phase = per[j]*range(0,stop=1,length=nphase) 
    lprob_phase = zeros(nphase)
    lprob_per[j] = -1e100
    # Loop over planet 3 phases:
    for i=1:nphase 
     # per param_names: mass ratio,phase,ecosw,esinw
      param_tmp = [log10(mratio),phase[i],0.01,0.01] 
      param3 = [best_p2;param_tmp] #concatenate 2 planet model to 3 planet model params
      per_cur = per[j]
      param1 = param3 .+ 100.0
      niter=0
      while maximum(abs.(param1 .- param3)) > tol && niter < 20
        param1 = param3
        fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:10];10^params[11];per_cur;params[12:end]],jmax,true),tt0,tt,weight,param3)
        param3 = fit.param
        niter+=1
        # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,param3,tt,sigtt,true,per_cur))
      end
      ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[param3[1:10];10^param3[11];per_cur;param3[12:end]],jmax,true)
      lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
      if lprob_phase[i] > lprob_best 
      # Check that best fit for phase is better than global best minimum
        lprob_best = lprob_phase[i]
        perbest = [fit.param[1:10];10^fit.param[11];per_cur;fit.param[12:end]]
      end
      if lprob_phase[i] > lprob_per[j] 
      # Check best fit over planet phases for this particular period
        lprob_per[j] = lprob_phase[i]
        param_per[1:nparam,j] = [fit.param[1:10];10^fit.param[11];per_cur;fit.param[12:end]]
      end
      # if j>1 && abs(lprob_p3[j] - lprob_p3[j-1])>5
      #   # Check that best fit for current period is close to that of previous period
      #   lprob_p3[j] = lprob_p3[j-1]
      #   param_p3[1:nparam,j] = [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
      # end
    end
    #println("Period: ",per[j]," log Prob: ",lprob_per[j]," Param: ",vec(param_per[1:nparam,j]))
  end
  println("Finished ",planet," planet fit w/ fixed period: ",perbest," in ",niter," iterations")
  
  writedlm(outfile,zip(per,lprob_per))
	
	best_per,err=global_fit(tt0,nplanet,ntrans,perbest,jmax,true)
  # @save outfile per lprob_per best_per lprob_best_per ntrans nplanet tt0 tt ttmodel sigtt
  pname=["mu_1","P_1","t01","ecos1","esin1",
          "mu_2","P_2","t02","ecos2","esin2",
          "mu_3","P_3","t03","ecos3","esin3"]
  #for i=1:nparam
  #  println(pname[i]," : ",best_per[i]," ± ",err[i])
  #end
  mean_mp=[best_per[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  mean_ecc=[sqrt(best_per[(iplanet-1)*5+4]^2 + best_per[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  ecc_errs=[sqrt(err[(iplanet-1)*5+4]^2 + err[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]

  results = string("results/mystery",planet,"_fitresults.txt")
  open(results,"w") do io
  	println(io,"Global Fit Results.",'\n',"per=[",per_in," - ",per_out,", length=",nper,"]")
  	for i=1:nparam
			println(io,pname[i],": ",best_per[i]," ± ",err[i])
  	end
    println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
    println(io,"Retrieved eccentricity:",'\n',mean_ecc,'\n'," ± ",ecc_errs)
  end
  fig=figure(figsize=(6,6))
  subplots_adjust(hspace=0.05,wspace=0.05)
  ax1=gca()
  lim=minimum(per),maximum(per)
  # xbin,xhist,xbin_square,hist_square=histogram(param,nbins)
  ax1.plot(per ./365.25,exp.(lprob_per .- maximum(lprob_per))) 
  ax1.axvline(per_guess/365.25,linestyle="--",color="black")
  ax1.text((per_guess/365.25) + .1,1.01,planet)
  ax1.set_xlabel("Planet Period Search Grid [years]")
  ax1.set_ylabel("Probability")
  ax1.minorticks_on()
  ax1.tick_params(which="both",direction="in")
  show()
  return #best_per
end

# function fit_planet3(filename::String,label::String,
#   jd1::Float64,jd2::Float64,jdsize::Int64,
#   perin::Float64,perout::Float64,nper::Int,nphase::Int,
#   sqrte::Bool=false,
#   addnoise::Bool=false,sigma::Float64=0.0,EM::Bool=true)

#   data1 = readdlm(filename)
#   nt1 = sum(data1[:,1] .== 1.0)
#   nt2 = sum(data1[:,1] .== 2.0)
#   tt1 = vec(data1[1:nt1,3])
#   tt2 = vec(data1[nt1+1:nt1+nt2,3])
  
#   if addnoise 
#     sigtt1 = data1[1:nt1,4]
#     sigtt2 = data1[nt1+1:nt1+nt2,4]
#   else
#     sigtt1 = ones(nt1)
#     sigtt2 = ones(nt2)
#   end

#   # Okay,let's do a linear fit to the transit times (third column):
#   function find_coeffs(tt,period,sigtt)
#     nt = length(tt)
#     x = zeros(2,nt)
#     x[1,1:nt] .= 1.0
#     x[2,1] = 0.0 
#     for i=2:nt
#       x[2,i] = x[2,i-1] + round((tt[i]-tt[i-1])/period) 
#     end
#     coeff,covcoeff = regress(x,tt,sigtt)
#     # println(tt,sigtt,std(sigtt))
#     return coeff,covcoeff
#   end

#   p1est = median(tt1[2:end] - tt1[1:end-1])
#   p2est = median(tt2[2:end] - tt2[1:end-1])

#   coeff1,covcoeff1 = find_coeffs(tt1,p1est,sigtt1)
#   coeff2,covcoeff2 = find_coeffs(tt2,p2est,sigtt2)

#   sigtt=[sigtt1;sigtt2] 
#   # @assert (sigtt[1] .* (24 * 3600) .= sigma)

#   t01 = coeff1[1]; per1 = coeff1[2]
#   t02 = coeff2[1]; per2 = coeff2[2]
#   t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) #best fit linear transit times w/o ttvs
#   t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
#   # Best-fit linear transit times:
#   tt0 = [t1;t2]
#   weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
#   # Actual transit times:
#   tt=[tt1;tt2]

#   # Okay,now let's do a 2-planet fit:
#   # param_names = mass ratio,period,initial transit time,e*cos(omega),e*sin(omega)
#   init_param = [3e-6,per1,t01,0.01,0.01,
#                 3e-6,per2,t02,0.01,0.01] 
#   println("Initial parameters: ",init_param)
#   #model = ttv_wrapper2(tt0,param)
#   # Set up data structure to hold planet properties,passed to TTVFaster
#   jmax = 5
#   data=init_param
#   p1=TTVFaster.Planet_plane_hk(data[1],data[2],data[3],data[4],data[ 5])
#   p2=TTVFaster.Planet_plane_hk(data[6],data[7],data[8],data[9],data[10])
#   time1 = collect(p1.trans0 .+ range(0,stop=nt1-1,length=nt1) .* p1.period)
#   time2 = collect(p2.trans0 .+ range(0,stop=nt2-1,length=nt2) .* p2.period)
#   # Initialize the computation of the Laplace coefficients:
#   ttv1 = zeros(nt1)
#   ttv2 = zeros(nt2)
#   # Need first call to TTVFaster,without optimizing
#   dummy=TTVFaster.compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2) 
#   # function plot_2planetfit(perin,perout,sigma)
#   #   clf()
#   #   scatter(time1,tt1.-t1)
#   #   plot(time1,ttv1)
#   #   scatter(time2,tt2.-t2,color="green")
#   #   plot(time2,ttv2)
#   #   name = string("IMAGES/2planetfitp",label,".png")
#   #   savefig(name)
#   # end

#   # Now,optimize 2-planet fit
#   per_cur = 11.86*365.25 #jupiter period in days,initial value
#   #res = optimize(chisquare2,param,method = :l_bfgs,iterations = 21)
#   ntrans = [nt1,nt2]
#   Nobs = sum(ntrans)
#   nplanet = 2
#   # create initial simplex? need function for this?
#   # result = optimize(f0,xcurr,NelderMead(initial_simplex=MySimplexer(),show_trace=true,iterations=1))
#   println("Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,sqrte))
#   # res = optimize(params -> chisquare(nplanet,ntrans,params,tt,sigtt),init_param) 
#   # init_param = res.minimizer
#   param1 = init_param .+ 100.0
#   while maximum(abs.(param1 .- init_param)) > 1e-5
#     param1 = init_param
#     res = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,sqrte),tt0,tt,weight,init_param)
#     init_param = res.param
#     # println("init_param: ",init_param)
#     # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt))
#   end
#   # res = optimize(params -> chisquare(tt0,nplanet,ntrans,params,tt,sigtt),init_param) 
#   # init_param = res.minimizer
#   # fit2 = curve_fit(ttv_wrapper2,tt0,tt,weight,param; show_trace=true)
#   println("Finished 2-planet fit: ",init_param)

#   # Now,let's add the 3rd planet:
#   ntrans = [nt1,nt2,2] #requires at least 2 transits for each planet (even if it doesnt transit)
#   nplanet = 3
#   #per = 11.86*365.25
#   # Grid of periods to search over:
#   per = 10 .^ range(log10(perin),stop=log10(perout),length=nper)
#   lprob_per = zeros(nper)
#   nparam = 15
#   param_per = zeros(nparam,nper)
#   lprob_best = -1e100 #global best fit
#   pbest = zeros(nparam)
#   # Shifting to simulated observation range to search over period grid
#   offset = (jd1 + jd2)/2 
#   for j=1:nper
#     phase = per[j]*range(0,stop=1,length=nphase) .+ offset 
#     lprob_phase = zeros(nphase)
#     lprob_per[j] = -1e100
#     for i=1:nphase #loops over jupiter phases
#       param_tmp = [1e-3,phase[i],0.01,0.01] # jupiter params: mass ratio,phase,ecosw,esinw
#       param3 = [init_param;param_tmp] #concatenate 2 planet model to 3 planet model params
#       per_cur = per[j] #sets jupiter period to global value
#       # fit = curve_fit(ttv_wrapper_fixper,tt0,tt,weight,param3) #optimizes fit w/ 3 planet model
#       # fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,true,per_cur),tt0,tt,weight,param3) 
#       # param3 = fit.param
#       param1 = param3 .+ 100.0
#       while maximum(abs.(param1 .- param3)) > 1e-5
#         param1 = param3
#         fit = curve_fit((tt0,param3) -> ttv_wrapper(tt0,nplanet,ntrans,[param3[1:11];per_cur;param3[12:end]],jmax,sqrte),tt0,tt,weight,param3)
#         param3 = fit.param
#         # println("init_param: ",param3)
#         # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,param3,tt,sigtt,true,per_cur))
#       end
#       ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[param3[1:11];per_cur;param3[12:end]],jmax,sqrte)
#       lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
#       if lprob_phase[i] > lprob_best # check that best fit for period is better than global best fit
#         lprob_best = lprob_phase[i]
#         pbest = [fit.param[1:11];per_cur;fit.param[12:14]]
#       end
#       if lprob_phase[i] > lprob_per[j] # checks best fit over all phases of jupiter for this particular period
#         lprob_per[j] = lprob_phase[i]
#         param_per[1:nparam,j] =  [fit.param[1:11];per_cur;fit.param[12:14]]
#       end
#     end
#     println("Period: ",per[j]," chi: ",lprob_per[j]," Param: ",vec(param_per[1:nparam,j]))
#   end
#   println("Finished 3-planet fit w/ fixed period: ",pbest)
  
#   # function plot_likelihood(perin,perout,sigma)
#   #   clf()
#   #   plot(per/365.25,exp.((lprob_per .-maximum(lprob_per)))) 
#   #   xlabel("Period of planet 3 [years]")
#   #   ylabel("Likelihood")
#   #   name = string("IMAGES/perlikelihood",label,".png")
#   #   savefig(name)
#   # end

#   #ttmodel=ttv_wrapper3(tt0,param3)
#   #res = optimize(chisquare3,param3,method = :l_bfgs,iterations = 21)
#   #  res = optimize(chisquare3,param3,method = :l_bfgs)
#   #  ttmodel=ttv_wrapper3(tt0,param3)

#   # fit = curve_fit(ttv_wrapper3,tt0,tt,weight,pbest)
#   fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,sqrte),tt0,tt,weight,pbest)
#   # ttmodel=ttv_wrapper3(tt0,pbest)
#   pbest_global = fit.param
#   ttmodel = ttv_wrapper(tt0,nplanet,ntrans,pbest_global,jmax,sqrte)
#   lprob_best= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
#   sigsys2 = 1e-6

#   println("Finished global 3-planet fit.")
#   println("Maximum: ",lprob_best," Param: ",pbest_global)

#   # function plot_3planetfit(perin,perout,sigma)
#   #   clf()
#   #   scatter(time1,tt1.-t1)
#   #   plot(time1,ttmodel[1:nt1].-t1)
#   #   scatter(time2,tt2.-t2,color="green")
#   #   plot(time2,ttmodel[nt1+1:nt1+nt2].-t2)
#   #   name = string("IMAGES/3planetfitp",label,".png")
#   #   savefig(name)
#   # end

#   pname = ["mu_1","P_1","t01","e1 cos(om1)","e1 sin(om1)",
#         "mu_2","P_2","t02","e2 cos(om2)","e2 sin(om2)",
#         "mu_3","P_3","t03","e3 cos(om3)","e3 sin(om3)"]

#   results = string("OUTPUTS/per_fit",label,"results.txt")
#   open(results,"w") do io
#     for i=1:nparam
#       println(io,pname[i],": ",pbest_global[i])
#     end
#   end
#   file = string("OUTPUTS/per_fit",label,"params.jld2")
#   @save file param_per lprob_per lprob_best pbest_global ntrans nplanet jd1 jd2 jdsize tt0 tt ttmodel sigtt perin perout nper nphase
#   # writedlm(results,pbest_global)
#     return lprob_best,pbest_global

#   # Now,search for Moon:
#   nparam = 18
#   deltaphi_cur = 2.312
#   deltaphi = range(dpin,stop=dpout,length=ndp)
#   lprob_dp = zeros(ndp)
#   param_dp = zeros(nparam,ndp)
#   lprob_best = -1e100 #global best fit
#   pbest_dp = zeros(nparam)
#   for j=1:ndp
#     lprob_dp[j] = -1e100 
#     param_tmp = [0.01,0.01,deltaphi[j]] # lunar params: t_s ,t_c ,deltaphi 
#     param4 = [pbest_per;param_tmp]
#     deltaphi_cur = deltaphi[j]
#     param1 = param4 .+ 100.0
#     while maximum(abs.(param1 .- param4)) > 1e-5
#       param1 = param4
#       fit = curve_fit((tt0,param4) -> ttv_wrapper(tt0,nplanet,ntrans,param4,jmax,false),tt0,tt,weight,param4)
#       param4 = fit.param 
#     end
#     ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[fit.param[1:17];deltaphi_cur],jmax,false)
#     lprob_dp[j]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
#     if lprob_dp[j] > lprob_best 
#       lprob_best = lprob_dp[j]
#       pbest_dp = [fit.param[1:17];deltaphi_cur]
#     end
#     # end
#     param_dp[1:nparam,j] = [fit.param[1:17];deltaphi_cur]
#     println("deltaphi: ",deltaphi[j]," chi: ",lprob_dp[j]," Param: ",vec(param_dp[1:nparam,j]))
#   end

#   fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,false),tt0,tt,weight,pbest_dp)
#   pbest_global = fit.param
#   ttmodel = ttv_wrapper(tt0,nplanet,ntrans,pbest_global,jmax,false)
#   lprob_best = (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
#   println("Finished lunar fit: ",lprob_best," ",pbest_global)

#   pname = ["mu_1","P_1","t01","e1 cos(om1)","e1 sin(om1)",
#             "mu_2","P_2","t02","e2 cos(om2)","e2 sin(om2)",
#             "mu_3","P_3","t03","e3 cos(om3)","e3 sin(om3)",
#             "tmax sin(phi0)","tmax cos(phi0)","deltaphi"]

#   results = string("OUTPUTS/moon_fit",label,"results.txt")
#   open(results,"w") do io
#     for i=1:nparam
#       println(io,pname[i],": ",pbest_global[i])
#     end
#   end
#   file = string("OUTPUTS/moon_fit",label,"params.jld2")
#   @save file pbest_per pbest_dp lprob_per lprob_dp lprob_best pbest_global ntrans nplanet jd1 jd2 jdsize tt0 tt ttmodel sigtt perin perout nper nphase dpin dpout ndp
#   # results = string("OUTPUTS/per_fit",label,"results.txt")
#   # #writedlm(results,pbest)
#   return lprob_best,pbest_global
# end

