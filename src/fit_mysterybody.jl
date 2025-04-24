# Julia v1.1
# if !@isdefined(TTVFaster)
#     include("TTVFaster/src/TTVFaster.jl")   
# end
# using Main.TTVFaster
# import Main.TTVFaster.ttv_wrapper
# import Main.TTVFaster.chisquare
# include("sim_times.jl") 
include("misc.jl")
# include("plot_likelihood.jl")
using TTVFaster,DelimitedFiles,JLD2,LsqFit,Statistics,Profile,PyPlot,CSV,DataFrames
rc("lines",linewidth=2)

# using PyPlot,Unitful,UnitfulAstro,LinearAlgebra
sigma=60; nyear=30#25#20#15
jd1=2.4332825e6
tref = 2430000
tol = 1e-5
per_guess = 11.88*365.25
per_in=1.6*365.25
per_out=22.0*365.25
nper, nphase = 200, 36
nmu=10
jmax = 5
planet="p3";obs= "fromEMB"
#sigma,nyear=parse(Int64,ARGS[1]),parse(Int64,ARGS[2])
datafile=string("INPUTS/EMBtt_",sigma,"s",nyear,"yrs.txt")
save_as_jld2=true
is_txt_file=true
# p4=jldopen("FITS/fromEMB/p4_fit30s30yrs.jld2","r")
#body,trans,tt,sigtt,tt0,noise=sim_obs_and_find_times(jd1,sigma,nyear,obs)
  function global_fit(tt,tt0,sigtt,nplanet,ntrans,init_param,jmax,EM) 
    Nobs=sum(ntrans)
    weight = ones(Nobs)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
  	  # Perform fit with best params, and calculate covariances for parameters: 
    fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM),tt0,tt,weight,init_param)
    covar=estimate_covar(fit)
    
    best_global = fit.param ;   nparam=length(best_global)
    err=[sqrt(covar[i,j]) for i=1:nparam, j=1:nparam if i==j ]
    ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_global,jmax,EM)
    lprob_best_global= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
    println("Finished global fit.")
    # println("New chi-square: ",chisquare(tt0,nplanet,ntrans,best_global,tt,sigtt,jmax,EM))
    # println("Maximum: ",lprob_best_global,'\n'," Param: ",best_global)
    return best_global,err
   end

function fit_mysteryplanet(;is_txt_file=true,is_jld_file=false) #fit_mysteryplanet(datafile::String,jd1::Float64,tref::Real,tol::Real,obs::String)#::Int,mratio::Float64,per_guess::Float64,per_in::Float64,per_out::Float64,nper::Int,nphase::Int,
  #outfile = string("FITS/",obs,"/mystery_",planet,"_fit",jd1,"JED.jld2")
  if is_txt_file
  outfile=string("FITS/mysteryplanet_",planet,sigma,"s",nyear,"yrs.jld2")
  #outfile = string("FITS/mystery_",planet,"_fit",jd1,".jld2")
  @assert isfile(datafile)
  println(datafile," loaded.")
  data1 = readdlm(datafile,Float64)
  tt = data1[:,3] .- tref
  sigtt = data1[:,4]
  end
  #if is_jld_file
   # ntrans=[49,31]
    #tt=p4["tt"]
    #sigtt=p4["sigtt"]
 # end
  nplanet_cond = 2
  nparam=10
  #ntrans=[sum(body .== i) for i=1:nplanet_cond]
  nt1 = sum(data1[:,1] .== 1.0)
  nt2 = sum(data1[:,1] .== 2.0)

#@show ntrans[1],ntrans[2],length(tt)
  # Actual transit times:
  ntrans=[nt1,nt2]
    Nobs = sum(ntrans)
 # nt1,nt2 = ntrans[1],ntrans[2]
  tt1,sigtt1 = tt[1:nt1], sigtt[1:nt1]
  tt2,sigtt2 = tt[nt1+1:nt1+nt2], sigtt[nt1+1:nt1+nt2]
  weight = ones(nt1+nt2)./ sigtt.^2 #assigns each data point stat weight d.t. noise = 1/σ^2
  # Estimate the periods of each conditioned planet:
  p1est = median(tt1[2:end] - tt1[1:end-1])
  p2est = median(tt2[2:end] - tt2[1:end-1])
  @show p1est,p2est
  # Okay,let's do a linear fit to the transit times:
  x1,t01,per1 = linear_fit(tt1,p1est,sigtt1)
  x2,t02,per2 = linear_fit(tt2,p2est,sigtt2)
  # Best-fit linear transit times without TTVs:
  # t01 = coeff1[1]; per1 = coeff1[2]
  # t02 = coeff2[1]; per2 = coeff2[2]
  t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) 
  t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
  tt0 = [t1;t2]

  subplot(111)
  plot(t1,tt1.-t1,".",color="green")
  plot(t1,tt1.-t1,color="black")
  plot(t2,tt2.-t2,".",color="salmon")
  plot(t2,tt2.-t2,color="black")
  clf()
  # Okay,now let's do a 2-planet fit:
    # param_names = mass ratio,period,initial transit time,e*cos(omega),e*sin(omega)
  init_param_guess = [3e-6,per1,t01,0.01,0.01,
                      3e-6,per2,t02,0.01,0.01] 
    # println("Initial parameters: ",init_param)

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
    @show init_param
    # Now,optimize 2-planet fit
    println("Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,true))
    param1 = init_param .+ 100.0
    niter = 0
    while maximum(abs.(param1 .- init_param)) > tol && niter < 20
      param1 = init_param
      res = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,init_param)
      init_param = res.param
      niter += 1
      println("init_param: ",init_param)
      println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,init_param,tt,sigtt,jmax,true))
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

		#best_p2,err=global_fit(tt,tt0,sigtt,nplanet,ntrans,init_param_guess,jmax,true)
    return best_p2,err
  end
  @time best_p2,err = fit_cond_planets(tt0,tt,sigtt,jmax,nplanet_cond,init_param_guess)
  # Now,let's add the 3rd planet:
  ntrans = [ntrans;2] #requires at least 2 transits for each planet (even if it doesnt transit)
  nplanet = nplanet_cond + 1
  nparam = length(best_p2) + 5
  # Grid of periods to search over:
  per = 10 .^ range(log10(per_in),stop=log10(per_out),length=nper)
  # want a grid of masses instead of assuming its value
  mu= range(log10(1e-8), stop=log10(1e-2),length=nmu)
  lprob_best = -1e100 #global best fit
  chisq=zeros(nparam)

  lprob_per = zeros(nper,nmu)#zeros(nper)
  perbest = zeros(nparam)
  per_cur = per_guess 
  param_per = zeros(nparam,nper,nmu)#zeros(nparam,nper)
  niter = 0
  # p3best = zeros(nparam)
  # lprob_p3 = zeros(np3,length(mu3))
  # param_p3=zeros(nparam,np3,length(mu3))
  # Loop over planet 3 masses
  for k=1:length(mu)
      mu_cur = mu[k]
    for j=1:nper
      phase = per[j]*range(0,stop=1,length=nphase) 
      lprob_phase = zeros(nphase)
      lprob_per[j,k] = -1e100
      # Loop over planet 3 phases:
      for i=1:nphase 
       # per param_names: mass ratio,phase,ecosw,esinw
        param_tmp = [mu_cur,phase[i],0.01,0.01] 
        param3 = [best_p2;param_tmp] #concatenate 2 planet model to 3 planet model params
        per_cur = per[j]
        param1 = param3 .+ 100.0
        niter=0
        while maximum(abs.(param1 .- param3)) > tol && niter < 20
          param1 = param3
          fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,[params[1:10];10^mu_cur;per_cur;params[12:end]],jmax,true),tt0,tt,weight,param3)
          param3 = fit.param
          niter+=1
          # println("New Initial chi-square: ",chisquare(tt0,nplanet,ntrans,param3,tt,sigtt,true,per_cur))
        end
        # @show fit.param
        ttmodel = ttv_wrapper(tt0,nplanet,ntrans,[param3[1:10];10^mu_cur;per_cur;param3[12:end]],jmax,true)
        lprob_phase[i]= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
        if lprob_phase[i] > lprob_best 
        # Check that best fit for phase is better than global best minimum
          lprob_best = lprob_phase[i]
          perbest = [fit.param[1:10];10^mu_cur;per_cur;fit.param[12:end]]
        end
        if lprob_phase[i] > lprob_per[j,k] 
        # Check best fit over planet phases for this particular period and mass
          lprob_per[j,k] = lprob_phase[i]
          param_per[1:nparam,j,k] = [fit.param[1:10];10^mu_cur;per_cur;fit.param[12:end]]
        end
        # if j>1 && abs(lprob_p3[j] - lprob_p3[j-1])>5
        #   # Check that best fit for current period is close to that of previous period
        #   lprob_p3[j] = lprob_p3[j-1]
        #   param_p3[1:nparam,j] = [fit.param[1:10];10^fit.param[11];p3_cur;fit.param[12:end]]
        # end
      end # phase loop
      # println("Period: ",per[j]," Mass-ratio: ",10^mu3[k]," log Prob: ",lprob_per[j,k])#" Param: ",vec(param_per[1:nparam,j]))
    end # per loop
  end # mass loop
  # println("Finished ",planet," planet fit w/ fixed period: ",perbest," in ",niter," iterations")
  #best_per,err=global_fit(tt,tt0,sigtt,nplanet,ntrans,perbest,jmax,true)  
    fit = curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,perbest)
  cov=estimate_covar(fit)
  err=[sqrt(cov[i,j]) for i=1:nparam, j=1:nparam if i==j ]
  best_per = fit.param
  ttmodel = ttv_wrapper(tt0,nplanet,ntrans,best_per,jmax,true)
  lprob_best_per= (1 - Nobs/2) * log(sum((tt-ttmodel).^2 ./sigtt.^2))
  # println("Finished global 3-planet fit.")
  println("New 3-planet chi-square: ",chisquare(tt0,nplanet,ntrans,best_per,tt,sigtt,jmax,true))
  println("Maximum: ",lprob_best_per," Param: ",best_per)
# writedlm(outfile,zip(per,lprob_per))
  # df=DataFrame(mu_1=param_per[1,:],P_1=param_per[2,:],t01=param_per[3,:],ecos1=param_per[4,:],esin1=param_per[5,:],
  #             mu_2=param_per[6,:],P_2=param_per[7,:],t02=param_per[8,:],ecos2=param_per[9,:],esin2=param_per[10,:],
  #             mu_3=param_per[11,:],P_3=param_per[12,:],t03=param_per[13,:],ecos3=param_per[14,:],esin3=param_per[15,:],
  #             lprob=lprob_per[:,i])
# CSV.write(grid,df)
	
# end
  # @save outfile per lprob_per best_per lprob_best_per ntrans nplanet tt0 tt ttmodel sigtt
  # pname=["mu_1","P_1","t01","ecos1","esin1",
  #         "mu_2","P_2","t02","ecos2","esin2",
  #         "mu_3","P_3","t03","ecos3","esin3"]
  #for i=1:nparam
  #  println(pname[i]," : ",best_per[i]," ± ",err[i])
  #end
  # mean_mp=[best_per[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  # mp_errs=[err[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  # mean_ecc=[sqrt(best_per[(iplanet-1)*5+4]^2 + best_per[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]
  # ecc_errs=[sqrt(err[(iplanet-1)*5+4]^2 + err[(iplanet-1)*5+4]^2) for iplanet=1:nplanet]

  # results = string("results/mystery",planet,"_fitresults.txt")
  # open(results,"w") do io
  # 	println(io,"Global Fit Results.",'\n',"per=[",per_in," - ",per_out,", length=",nper,"]")
  # 	for i=1:nparam
	# 		println(io,pname[i],": ",best_per[i]," ± ",err[i])
  # 	end
  #   println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
  #   println(io,"Retrieved eccentricity:",'\n',mean_ecc,'\n'," ± ",ecc_errs)
  # end

  if save_as_jld2
  @save outfile per lprob_per best_per lprob_best_per ntrans nplanet tt0 tt ttmodel sigtt nphase param_per
  end
  return best_per
end
#best_per = fit_mysteryplanet()
cmap=plt.cm.get_cmap("plasma")
col_cycler=plt.cycler("color",cmap(range(0,1,length=nmu)))
ls_cycler=plt.cycler("linestyle",repeat(["--",":"],5))
wp3=jldopen(string("FITS/mysteryplanet_",planet,sigma,"s",nyear,"yrs.jld2"),"r")

# wp3=jldopen("FITS/fromEMB/widep3_fit100s30yrs.jld2","r");
# sigma=100
data=wp3["tt"]
label_cycle = plt.cycler(label=["set {n}" for n in 1:4])
# println(xprob(wp3["lprob_p3"],data))
# logL=xprob(actual_logL(wp3["lprob_p3"][:,1],data))
# print("lnL",logL)
mu= range(log10(1e-8), stop=log10(1e-2),length=nmu)
#function plot_mu()# plot grid 

   fig,ax=subplots(3,2,figsize=(10,4.5))
  ax1=plt.subplot2grid((3,2),(0,0),rowspan=3)
# ax1.text(.5,.5,"[0,0]")
  ax2=plt.subplot2grid((3,2),(0,1))
  # ax2.text(.5,.5,"[0,2]")
  ax3=plt.subplot2grid((3,2),(1,1))
  ax4=plt.subplot2grid((3,2),(2,1))
  # ax3.text(.5,.5,"[1,2]")
  # show()
  # ax4=plt.subplot2grid((3,2),(2,1))
  # ax5=plt.subplot2grid((3,3),(0,2))
  # ax6=plt.subplot2grid((3,3),(1,2))
  ax1.set_prop_cycle(col_cycler+ls_cycler)
  ax2.set_prop_cycle(col_cycler+ls_cycler)
  ax3.set_prop_cycle(col_cycler+ls_cycler)
  ax1.axvline(11.86,linestyle="-",color="black")
  ax1.axvline(1.88,linestyle="-",color="black")
  ax2.axvline(11.86,linestyle="-",color="black")

  ax2.axvline(1.88,linestyle="-",color="black")
  ax3.axvline(11.86,linestyle="-",color="black")
  ax4.set_frame_on(false)
  ax4.set_xticks([])
  ax4.set_yticks([])
inds = argmax(wp3["lprob_per"])
@show inds
  # ax3=fig.add_axes([0.78,0.2,0.2,0.3])
  # ax3.set_prop_cycle(col_cycler+ls_cycler)
  handles=[]
  for i=1:nmu
    ax2.plot(wp3["per"]./365.35,xprob(actual_logL(wp3["lprob_per"][:,i],data)))
    ax3.plot(wp3["per"]./365.35,xprob(wp3["lprob_per"][:,i],data,wp3["lprob_per"][inds[1],inds[2]]))
   h = ax1.plot(wp3["per"]./365.35,actual_logL(wp3["lprob_per"][:,i],data),label=string(round(mu[i],sigdigits=2)))
    push!(handles,h)
  end
  # ax3.set_xlim(10,13)
  # ax3.tick_params(left="false",labelleft="false",right="true",labelright="true")
  # fig.legend(labels=string.(round.(collect(mu),sigdigits=2)),handles=handles)
  fig.legend(loc="upper right",title=string(L"$\gamma$ values"),fontsize="medium",title_fontsize="large",bbox_to_anchor=(0.09,0.25,0.85,.102),ncol=2)
  # fig.legend(loc="upper right",title=string(L"$\gamma$ values"),fontsize="medium",title_fontsize="large",ncol=4)
  fac(true_per)=true_per + true_per/100
  ax1.text(fac(1.88),-170,"Mars",color="black")
  ax1.text(fac(11.88),-170,"Jupiter",color="black")
  ax4.text(0.05,0.2,string(L"$σ_{obs}$","=",sigma," s",'\n',L"$n_{years}$","=",nyear,'\n',L"$m_p/ M_{\odot}$=",L"$10^γ$"),color="black",fontsize="large")
  # ax[1].text(3,10000,string(L"$m_p/ M_{\odot}$=",L"$10^γ$"),color="black")
  ax1.set_ylim(-220,-160)
  ax1.minorticks_on()
  ax2.minorticks_on()
  ax3.minorticks_on()

  ax3.set_ylabel("Relative Prob.",fontsize="medium")
  ax2.set_ylabel("Relative Prob.",fontsize="medium")
  ax3.set_title("With respect to all mass-ratios")
  ax1.set_ylabel(string("ln",L"$\mathcal{L}$"),fontsize="x-large")
  ax3.tick_params(top=true,direction="in")
  ax1.set_xlabel("Orbital Period [years]",fontsize="x-large")
  ax3.set_xlabel("Orbital Period [years]",fontsize="medium")
 # fig.supxlabel("Orbital Period [years]",fontsize="x-large")
  # fig.suptitle("Difference between actual and approximate logL in blind search")
  tight_layout()
  savefig("2025/2025actual_logL_zoom_jup_norm_60_30_.png",dpi=150)
#end
#plot_mu()
#show()
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

#=
  fig=figure(figsize=(6,6))
  subplots_adjust(hspace=0.05,wspace=0.05)
  ax1=gca()
  lim=minimum(per),maximum(per)
  # xbin,xhist,xbin_square,hist_square=histogram(param,nbins)
  ax1.plot(per ./365.25,exp.(lprob_per .- maximum(lprob_per))) 
  ax1.axvline(per_guess/365.25,linestyle="--",color="black")
  ax1.text((per_guess/365.25) + .1,1.01,planet)
  ax1.set_xlabel("Planet Period Search Grid [years]")
  ax1.set_ylabel("Relative Probability")
  ax1.minorticks_on()
  ax1.tick_params(which="both",direction="in")
  show()
 =#
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

