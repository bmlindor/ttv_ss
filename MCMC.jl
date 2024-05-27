include("bounds.jl")
#include("CGS.jl")
include("misc.jl")
using TTVFaster,DelimitedFiles,JLD2,LaTeXStrings,PyPlot
using Statistics,StatsBase,MCMCDiagnostics

# Run a Markov chain:
function MCMC(foutput::String,param::Array{Float64,1},lprob_best::Float64,nsteps::Int64,nwalkers::Int64,nplanet::Int64,ntrans::Array{Int64,1},tt0::Array{Float64,1},tt::Array{Float64,1},sigtt::Array{Float64,1},use_sigsys::Bool,EM::Bool) 
  nparam = length(param)
  println(nparam," parameters from fit: ",param)
  println("Maximum log Prob from fit: ",lprob_best)
	println("No Walkers: ",nwalkers," No. steps: ",nsteps)
  jmax = 5
  errors = [1e-7,1e-5,1e-5,1e-2,1e-2,
            1e-7,1e-5,1e-5,1e-2,1e-2,
            1e-6,1e-1,1e-1,1e-2,1e-2]
  pname = ["mu_1","P_1","t01","e1 cos(om1)","e1 sin(om1)",
          "mu_2","P_2","t02","e2 cos(om2)","e2 sin(om2)",
          "mu_3","P_3","t03","e3 cos(om3)","e3 sin(om3)"]
  if nplanet==3
    errors = errors
    pname = pname
  elseif nplanet==4
    errors = [errors[1:10];1e-7;1e-5;1e-5;1e-2;1e-2;errors[11:end]]
    pname = [pname[1:end];"mu_4";"P_4";"t04";"e4 cos(om4)";"e4 sin(om4)"]
  elseif nplanet==5
    errors=[errors[1:10];1e-7;1e-1;1e-1;1e-2;1e-2;errors[11:end];1e-6;1e-1;1e-1;1e-2;1e-2]
    pname = [pname[1:end];"mu_4";"P_4";"t04";"e4 cos(om4)";"e4 sin(om4)";"mu_5";"P_5";"t05";"e5 cos(om5)";"e5 sin(om5)"]
  elseif nplanet==2
    errors=errors[1:10]
    pname=pname[1:10]
  end
  if mod(nparam,5) == 3 
    moon_errors = [1e-4,1e-4,1e-5]
    moon_name = ["tcosϕ","tsinϕ","Δϕ"]
    append!(errors,moon_errors)
    append!(pname,moon_name)
  end
  # Initialize walkers:
  if use_sigsys
    pname=[pname;"σ_sys"]
    par_trial = [param[1:end];1e-8]
    nsize = nparam+1
  else
    par_trial = param[1:end]
    nsize = nparam
  end
  # Number of walkers should be at least 3 times number of parameters:
  @assert (nwalkers >= 3*nparam)  
  # Set up arrays to hold the results:
  par_mcmc = zeros(nwalkers,nsteps,nsize)
  lprob_mcmc = zeros(nwalkers,nsteps)
	thinning=1000
	thin_chain = zeros(nwalkers, div(nsteps,thinning), nsize)
	thin_lprob = zeros(nwalkers,div(nsteps,thinning))

  function calc_lprior(param) 
  # We will place a joint prior on eccentricity vector
  # such that each planet has an eccentricity which lies between
  # 0 and emax2,with a gradual decrease from emax1 to emax2:
  # Eccentricity priors:
    emax1 = 0.2; emax2 = 0.3
    # Define -log(prior):
    lprior = 0.0
    # Loop over planets:
    for iplanet=1:nplanet
    # Place prior on eccentricities:
      ecc = sqrt(param[(iplanet-1)*5+4]^2+param[(iplanet-1)*5+5]^2)
      lprior_tmp,dpdx = log_bounds_upper(ecc,emax1,emax2)
      lprior += lprior_tmp
      lprior += -log(ecc) 
    end
    for iplanet=1:nplanet-1
    # The periods of the planets should be ordered from least to greatest:
        if param[(iplanet-1)*5+2] > param[iplanet*5+2]
          lprior += -Inf
        end
    end
    for iplanet=1:nplanet
    # The masses should be positive:
      if param[(iplanet-1)*5+1] < 0 #unsure if this is correct implement
        lprior += -Inf
      end
    end
		if mod(nparam,5) == 3
    # Force the deltaphi to be between 0 and pi (to account for aliasing):
      dpmin = 0.0; dpmax = pi
      deltaphi = param[18]
      while deltaphi < dpmin
        deltaphi += 2pi
      end
      while deltaphi > 2pi
        deltaphi -= 2pi 
      end
      if deltaphi > dpmax
        deltaphi = 2pi - deltaphi
      end
      param[18] = deltaphi
      # if deltaphi < dpmin || deltaphi > dpmax
      #   lprior += -Inf
      # end
    end
    if use_sigsys
    # The systematic uncertainty should be positive:
      sigsys = param[end]
      if sigsys < 0 
        lprior += -Inf
      end
    end
    return lprior
  end

  Nobs = sum(ntrans) - 2*(nplanet - 2)
  for j=1:nwalkers
  # Select from within uncertainties:
    lprob_trial = -1e100
  # Only initiate models with reasonable Log Prob values:
    while lprob_trial < lprob_best - 1000 #since logProb, maybe 1000 is too large
      par_trial[1:nparam] .= param .+ errors .* randn(nparam) 
      if use_sigsys
        par_trial[nparam+1] = 1e-8 .* abs(randn())
      end
      lprob_trial = calc_lprior(par_trial)
      # println("Calculated Log Prior: ",lprob_trial)
      if lprob_trial > -Inf
        model = ttv_wrapper(tt0,nplanet,ntrans,par_trial[1:nparam],jmax,EM)
        # println(length(tt)," ",length(model))
        if use_sigsys
    # ll = (-0.5 * sum((tt-model).^2 ./(sigtt.^2 .+ sigsys.^2) .+ log.(sigtt.^2 .+ sigsys.^2)))
         lprob_trial += (-0.5 * sum((tt-model).^2 ./(sigtt.^2 .+ par_trial[end].^2) .+ log.(sigtt.^2 .+ par_trial[end].^2)))
        else
    # ll = log(sum((tt-model).^2 ./sigtt.^2))
          lprob_trial += log(sum((tt-model).^2 ./sigtt.^2))*(1 - Nobs/2) #mostly useful for grid search
        end
      end
      # println("Trial Log Prob: ",lprob_trial)
    end
    lprob_mcmc[j,1]=lprob_trial
    par_mcmc[j,1,:]=par_trial
    # println("Success: ",par_trial,lprob_trial)
  end

  # Initialize scale length & acceptance counter:
  ascale = 2.0
  accept = 0
  # Next,loop over steps in markov chain:
  for i=2:nsteps
    for j=1:nwalkers
      ipartner = j
  # Choose another walker to 'breed' a trial step with:
      while ipartner == j
        ipartner = ceil(Int,rand()*nwalkers)
      end
  # Now,choose a new set of parameters for the trial step:
      z=(rand()*(sqrt(ascale)-1.0/sqrt(ascale))+1.0/sqrt(ascale))^2 # draws from 1=sqrt(z....)
      par_trial=vec(z*par_mcmc[j,i-1,:]+(1.0-z)*par_mcmc[ipartner,i-1,:])
  # Compute model & Log Prob for trial step:  
      lprob_trial = calc_lprior(par_trial)  
      if lprob_trial > -Inf
          model_trial = ttv_wrapper(tt0,nplanet,ntrans,par_trial[1:nparam],jmax,EM)
    # ll = -.5 *sum((tt-model_trial).^2 ./(sigtt.^2 .+ par_trial[end]) .+ log.(sigtt.^2 .+ par_trial[end]))
    # ll =  log(sum((tt-model_trial).^2 ./sigtt.^2))
    # lprob_trial = calc_lprior(par_trial) + ll*(1 - Nobs/2) 
        if use_sigsys
         lprob_trial += (-0.5 * sum((tt-model_trial).^2 ./(sigtt.^2 .+ par_trial[end].^2) .+ log.(sigtt.^2 .+ par_trial[end].^2)))
        else
          lprob_trial += log(sum((tt-model_trial).^2 ./sigtt.^2))*(1 - Nobs/2) 
        end 
      end   
  # Next,determine whether to accept this trial step:
	# Calculate ratio between Log Prob of trial step and previous step:
      alp = z^(nsize-1)*exp((lprob_trial - lprob_mcmc[j,i-1]))
			#logratio= log(z)*(nsize-1)+lprob_trial - lprob_mcmc[j,i-1]
			#if log(rand()) < logratio # if rejected, last values equal new values
			#  par_mcmc[j,i,:] = par_trial
      #  lprob_mcmc[j,i] = lprob_trial
			#end
     # if rand() < 0.0001 
     #   println("Step: ",i," Walker: ",j," Trial Log Prob: " ,lprob_trial," Prob: ",alp," Frac: ",accept/(mod(i-1,1000)*nwalkers+j))
     # end
      if alp >= rand()
  # If step is accepted,add it to the chains!
        par_mcmc[j,i,:] = par_trial
        lprob_mcmc[j,i] = lprob_trial
        accept = accept + 1
      else
  # If step is rejected,then copy last step:
        par_mcmc[j,i,:] = par_mcmc[j,i-1,:]
        lprob_mcmc[j,i] = lprob_mcmc[j,i-1]
      end
    end
    if mod(i,thinning) == 0
      println("Number of steps: ",i," Acceptance Rate: ",accept/(thinning*nwalkers))
			#thin_chain[j,div(i,thinning),:] = par_mcmc[j,i,:]
      #thin_lprob[j,div(i,thinning)] = lprob_mcmc[j,i]
 		end
  end

  # Now,determine time of burn-in by calculating first time median is crossed:
  iburn = 0
  for i=1:nsize
    med_param=median(par_mcmc[1:nwalkers,1:nsteps,i])
    for j=1:nwalkers
      istep=2
      while (par_mcmc[j,istep,i] > med_param) == (par_mcmc[j,istep-1,i] > med_param) && (istep < nsteps)
        istep=istep+1
      end
      if istep >= iburn
        iburn = istep
      end
    end
  end
  println("Burn-in Number (ends): ",iburn)

  # Calculate the minimum number of effective samples over all parameters:
  samplesize = zeros(nparam)
  for j=1:nwalkers
    for i=1:nparam
      samplesize[i] += effective_sample_size(par_mcmc[j,:,i])
    end
  end
  indepsamples = minimum(samplesize)
  # println("Independent Sample Size: ",indepsamples)

	# Calculate the auto-correlation function for each parameter:
	#acf=zeros(nparam)
	#for i=1:nparam
		#acf=autocorrelation(par_mcmc[:,iburn:nsteps,i],, )
	#end
	#println("Auto-Cor Function: ",acf)

  # Find mean and standard deviation of posteriors after burn-in:
	mean_posteriors=[mean(par_mcmc[:,iburn:nsteps,i]) for i=1:nparam]
	# println("Mean posteriors after burn-in: ",mean_posteriors)
  sigtot=[sqrt((mean(vec(par_mcmc[:,iburn:nsteps,end])).*3600*24)^2 + (mean(sigtt).*3600*24)^2) ]
	# println("Total timing uncertainty of: ",sigtot," secs.")
  for i=1:nparam
   avg[i],minus1sig[i],plus1sig[i]=quantile(vec(par_mcmc[:,iburn:end,i]),[0.5,0.1587,0.8413])
   println(pname[i]," = ",avg[i]," + ",abs(plus1sig[i]-avg[i])," _ ",abs(avg[i]-minus1sig[i]))
  end
# Find percentage of walkers where diff. between median and quantile value is >100
  # bad_walk=[]
  # for i in 1:nwalkers
  #   for j in 1:nparam
  #     walker_med,walker_quant=quantile!(par_mcmc[i,jldmc["iburn"]+1:end,j],[0.5,0.9])
  #     walk_start=par_mcmc[i,jldmc["iburn"]+1,j] 
  #     walk_end = par_mcmc[i,jldmc["iburn"]+1,j]
  #     ratio = walk_end/walk_start
  #     walker_prob=median(lprob_mcmc[i,jldmc["iburn"]+1:end])
  #     if abs(walk_end-walk_start)/walk_start > 0.1
  #       #abs(walker_med-walker_end)>30
  #       # println(i," ",walker_prob[i])
  #       append!(bad_walk,i)
  #     end
  #   end
  #     # If systematic uncertainty > injected uncertainty, reject
  #   # if median(par_mcmc[i,jldmc["iburn"]:end,end]).*3600*24 >= sigma
  #   #   println("Reject results?")
  #   #   append(bad_walk,i)
  #   # end
  # end
  # println("Bad walkers: ",bad_walk)

  # Plot traces
  # for i=2:nparam
  #   for j=1:i-1
  #     scatter(vec(par_mcmc[1:nwalkers,iburn:end,i]),vec(par_mcmc[1:nwalkers,iburn:end,j]))
  #     xlabel(pname[i])
  #     ylabel(pname[j])
  #        println("Hit return to continue")
  #     read(STDIN,Char)
  #     clf()
  #   end
  # end
  
  #open(fresults,"w") do io
	#println(io,"MC chi-square after burn-in: ",chi2)
  #for i=1:nparam
  #  println(io,pname[i], mean(vec(par_mcmc[:,iburn:nsteps,i]))," ± ",std(vec(par_mcmc[:,iburn:nsteps,i])))
  #end
  #  println(io,"Retrieved Earth masses:",'\n',mean_mp,'\n'," ± ",mp_errs)
  #  println(io,"Retrieved eccentricity:",'\n',mean_ecc,'\n'," ± ",ecc_errs)
	#	println(io,"Retrieved σ_tot [secs] : ",sigtot)
  #end
  #println("Saved in ",foutput)
  mcmcfile = string(foutput)
  @save mcmcfile par_mcmc lprob_mcmc param nwalkers nsteps iburn indepsamples pname
  return lprob_mcmc #, param, nwalkers, nsteps, accept, iburn, indepsamples
end

# Retrieve MCMC results after burn-in, remove bad walkers
function mc_vals(sigma::Real,nyear::Real,grid_type_nplanet::String,case_num=Int,include_moon::Bool=false)
  EM=true
  if case_num==1 && isfile(string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile=string("FITS/fromEMB/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
  elseif case_num==2 && isfile(string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile=string("FITS/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
  else
      return println("MCMC file for case ",case_num," with ",grid_type_nplanet," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  if include_moon
    EM=false
  end
  # println("Added moon to model? ",include_moon)
  jldmc=jldopen(String(mcfile),"r")
  jldfit=jldopen(String(fitfile),"r")
  nwalkers,nsteps=jldmc["nwalkers"],jldmc["nsteps"]
  iburn,samples=jldmc["iburn"], jldmc["indepsamples"]
  par_mcmc=jldmc["par_mcmc"]; lprob_mcmc=jldmc["lprob_mcmc"]  ; param=jldmc["param"]
  pname=jldmc["pname"]
  tt0,tt,ttmodel,sigtt=jldfit["tt0"],jldfit["tt"],jldfit["ttmodel"],jldfit["sigtt"]
  nplanet,ntrans=jldfit["nplanet"],jldfit["ntrans"]
  nt1,nt2=jldfit["ntrans"][1],jldfit["ntrans"][2]
  jmax=5
  # weight=ones(nt1+nt2)./ sigtt.^2 
  nparam=length(pname)
  prob_max=maximum(lprob_mcmc[:,iburn:end])
  sigsys=round((median(vec(par_mcmc[:,iburn:end,end]))).* 3600*24,sigdigits=3)
  sigsys_err=(std(vec(par_mcmc[:,iburn:end,end]))).* 3600*24
  sigtot=round(sqrt(sigsys^2 + sigma^2),sigdigits=4)
  # function plot_trace()
  fig, axs = plt.subplots(4,nplanet,figsize=(3*nplanet,nplanet*3))
  figtitle=string("MC Traces for ",sigma," s;",nyear," yr simulations of Venus and EMB")
  fig.suptitle(figtitle)
  count=0
  for j=1:nwalkers
      axs[1,1].plot(par_mcmc[j,iburn:nsteps,end].*24*3600,lprob_mcmc[j,iburn:nsteps])
      axs[1,2].plot(lprob_mcmc[j,iburn:nsteps])
      # axs[1,2].plot(par_mcmc[j,iburn:nsteps,end].*24*3600)
      for iplanet=1:nplanet
          axs[2,iplanet].plot(par_mcmc[j,iburn:end,(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH)
          axs[3,iplanet].plot(par_mcmc[j,iburn:end,(iplanet-1)*5+2])
          axs[4,iplanet].plot(par_mcmc[j,iburn:end,(iplanet-1)*5+4])
  #         count+=1
      end
  end
  # axs[1,4].set_ylabel("log Prob")
  axs[1,1].set_xlabel(L"$σ_{sys}$ [s]")
  axs[1,1].set_ylabel("log Prob")
  axs[1,2].set_ylabel(L"$σ_{sys}$ [s]")
  for iplanet=1:nplanet
      axs[2,iplanet].set_ylabel(L"$M [M_{\oplus}]$")
      axs[3,iplanet].set_ylabel(pname[(iplanet-1)*5+2])
      axs[4,iplanet].set_ylabel(pname[(iplanet-1)*5+4])
      if iplanet==3 || iplanet==4
        plt.delaxes(ax=axs[1,iplanet])
      end
  end

  # tight_layout()
  title=string("IMAGES/trace/case",case_num,grid_type_nplanet,"-",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  # end
  vals=jldmc["par_mcmc"][:,jldmc["iburn"]:end,:]#,sigdigits=6)
  reduced_chisq, BIC,chisq=round.(calc_BIC(jldmc["lprob_mcmc"][:,jldmc["iburn"]:jldmc["nsteps"]],jldfit["tt0"],jldfit["tt"],jldfit["sigtt"],jldfit["nplanet"],jldfit["ntrans"],vals),sigdigits=6)
  # @show BIC
	avg=zeros(nparam)
  med=zeros(nparam)
  low=zeros(nparam)
  errors=zeros((2,nparam))
  high=zeros(nparam)
  # st_dev=zeros(nparam)

  for i=1:nparam
   med[i],low[i],high[i]=quantile(vec(par_mcmc[:,iburn:end,i]),[0.5,0.1587,0.8413])
   avg[i]=mean(vec(par_mcmc[:,iburn:end,i]))
   errors[1,i]=med[i]-low[i]; errors[2,i]=high[i]-med[i]
    # st_dev[i]=std(vec(par_mcmc[:,iburn:end,i]))
   # println(pname[i]," = ",avg[i]," + ",abs(plus1sig[i]-avg[i])," _ ",abs(avg[i]-minus1sig[i]))
  end
  masses=[med[i-4] for i in 1:length(param) if i%5==0] .*CGS.MSUN/CGS.MEARTH
  # ecc=[calc_ecc(med[i-1],med[i]) for i in 1:length(param) if i%5==0] 
	# periods=[med[i-3] for i in 1:length(param) if i%5==0]
	
  println("Retrieved values.")
  println("M_p[M⊕]= ",masses)#" + ",masses.-mass_high," - ",masses.-mass_low)
  # println("std(M_p)= ",mass_errs)
  # # println("Per [d]= ",periods)#," +/- ",per_errs)
  # println("eccen. =",ecc)#," +/- ",ecc_errs)
  # println("σsys[s]= ",sigsys," +/- ",sigsys_err)
  # println("σtot[s]= ",sigtot)
  return med,errors
end

function mc_table(sigma::Real,nyear::Real,options,include_moon::Bool=false)
  obs=options[1]; fit_type_nplanet=options[2]; #bestfit=options[3]
  EM=true
  if obs=="fromEMB"
    mcfile=string("MCMC/fromEMB/",fit_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    mcfile2=string("MCMC/fromEMB/p2_mcmc",sigma,"s",nyear,"yrs.jld2")
    mcfile3=string("MCMC/fromEMB/p3_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile=string("FITS/fromEMB/",fit_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    fitfile2=string("FITS/fromEMB/p2_fit",sigma,"s",nyear,"yrs.jld2")
    fitfile3=string("FITS/fromEMB/p3_fit",sigma,"s",nyear,"yrs.jld2")
    label="EMB";case=1
    low_lim=-6.5;high_lim=6.5
  elseif obs=="fromEV"
    mcfile=string("MCMC/",fit_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    mcfile2=string("MCMC/p2_mcmc",sigma,"s",nyear,"yrs.jld2")
    mcfile3=string("MCMC/p3_mcmc",sigma,"s",nyear,"yrs.jld2")
    mcfile4=string("MCMC/widep4_mcmc",sigma,"s",nyear,"yrs.jld2")
    f4=jldopen(String(fitfile4),"r")
    mc4=jldopen(String(mcfile4),"r")
    label="Earth";case=2
    low_lim=-9;high_lim=9
    avg4=[mean(vec(mc4["par_mcmc"][:,mc4["iburn"]:end,i])) for i=1:21]
    sigsys4=round(avg4[end].*24*3600,sigdigits=3)
  # end
  #   return  println("FITS file for ",sim," with ",fitmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
    end
    if include_moon
      EM=false
    end

  f=jldopen(String(fitfile),"r")
  f2=jldopen(String(fitfile2),"r")
  f3=jldopen(String(fitfile3),"r")
  mc=jldopen(String(mcfile),"r")
  mc2=jldopen(String(mcfile2),"r")
  mc3=jldopen(String(mcfile3),"r")
  tt,tt0,sigtt,ttmodel=f["tt"],f["tt0"],f["sigtt"],f["ttmodel"]
  pname=mc["pname"];nparam=length(pname)

  vals=mc["par_mcmc"][:,mc["iburn"]:end,:]#,sigdigits=6)
  vals2=mc2["par_mcmc"][:,mc2["iburn"]:end,:]#,sigdigits=6)
  vals3=mc3["par_mcmc"][:,mc3["iburn"]:end,:]#,sigdigits=6)

  avg=[quantile(vec(mc["par_mcmc"][:,mc["iburn"]:end,i]),0.5) for i=1:nparam]#,sigdigits=6)
  avg2=[quantile(vec(mc2["par_mcmc"][:,mc2["iburn"]:end,i]),0.5) for i=1:11]#,sigdigits=6)
  avg3=[quantile(vec(mc3["par_mcmc"][:,mc3["iburn"]:end,i]),0.5) for i=1:16]#,sigdigits=6)

  low=round.([quantile(vec(mc["par_mcmc"][:,mc["iburn"]:end,i]),0.1587) for i=1:nparam],sigdigits=6)
  low2=round.([quantile(vec(mc2["par_mcmc"][:,mc2["iburn"]:end,i]),0.1587) for i=1:11],sigdigits=6)
  low3=round.([quantile(vec(mc3["par_mcmc"][:,mc3["iburn"]:end,i]),0.1587) for i=1:16],sigdigits=6)

  high=round.([quantile(vec(mc["par_mcmc"][:,mc["iburn"]:end,i]),0.8413) for i=1:nparam],sigdigits=6)
  high2=round.([quantile(vec(mc2["par_mcmc"][:,mc2["iburn"]:end,i]),0.8413) for i=1:11],sigdigits=6)
  high3=round.([quantile(vec(mc3["par_mcmc"][:,mc3["iburn"]:end,i]),0.8413) for i=1:16],sigdigits=6)

  # prob=quantile(exp.(mc["lprob_mcmc"][mc["iburn"]:mc["nsteps"]]),0.5);#prob_max = maximum(exp.(mc["lprob_mcmc"][mc["iburn"]:mc["nsteps"]]))
  #  prob2=quantile(exp.(mc2["lprob_mcmc"][mc2["iburn"]:mc2["nsteps"]]),0.5);#prob_max2 = maximum(exp.(mc2["lprob_mcmc"][mc2["iburn"]:mc2["nsteps"]]))
  #   prob3=quantile(exp.(mc3["lprob_mcmc"][mc3["iburn"]:mc3["nsteps"]]),0.5);#prob_max3 = maximum(exp.(mc3["lprob_mcmc"][mc3["iburn"]:mc3["nsteps"]]))
  #println(" median Prob: ",prob,"      maximum Prob: ",prob_max)
  #chi2_avg = chi_mcmc(tt0,nplanet,ntrans,mean_posteriors,tt,sigtt,jmax,EM)
  reduced_chi,BIC,chi=round.(calc_BIC(mc["lprob_mcmc"][:,mc["iburn"]:mc["nsteps"]],f["nplanet"],f["ntrans"],vals),sigdigits=6)
  reduced_chi2,BIC2,chi2=round.(calc_BIC(mc2["lprob_mcmc"][:,mc2["iburn"]:mc2["nsteps"]],f2["nplanet"],f2["ntrans"],vals2),sigdigits=6)
  reduced_chi3,BIC3,chi3=round.(calc_BIC(mc3["lprob_mcmc"][:,mc3["iburn"]:mc3["nsteps"]],f3["nplanet"],f3["ntrans"],vals3),sigdigits=6)

  # scatter1=(ttvmodel1.-ttv1)
  # scatter2=(ttvmodel2.-ttv2)
  # println("Venus Peak amplitude of O-C: ", maximum(scatter1))
  # println(label," Peak amplitude of O-C: ", maximum(scatter2))
  # println("Mean: ",mean(scatter1)," ",mean(scatter2))
  # println("Std: ",std(scatter1)," ",std(scatter2))

  #sigsys=round(avg[end].*24*3600,sigdigits=3)
  #sigsys2=round(avg2[end].*24*3600,sigdigits=3)
  #sigsys3=round(avg3[end].*24*3600,sigdigits=3)
  parname=[L"$\mu_1 \times 10^{-6}$",L"$P_1$ [days]",L"$t_{0,1}$ [days]",L"$e_1 \cos{\omega_1}$",L"$e_1 \sin{\omega_1}$",
            L"$\mu_2 \times 10^{-6}$",L"$P_2$ [days]",L"$t_{0,2}$ [days]",L"$e_2 \cos{\omega_2}$",L"$e_2 \sin{\omega_2}$",
            L"$\mu_3$",               L"$P_3$ [days]",L"$t_{0,3}$ [days]",L"$e_3 \cos{\omega_3}$",L"$e_3 \sin{\omega_3}$",
            L"$\mu_4 \times 10^{-6}$",L"$P_4$ [days]",L"$t_{0,4}$ [days]",L"$e_4 \cos{\omega_4}$",L"$e_4 \sin{\omega_4}$",
            L"$\mu_5$",               L"$P_5$ [days]",L"$t_{0,5}$ [days]",L"$e_4 \cos{\omega_5}$",L"$e_5 \sin{\omega_5}$",
            L"$t_{max} \sin{\phi_0}$",L"$t_{max} \cos{\phi_0}$",L"$\Delta \phi$ [rad]",L"$\sigma_{sys}^2$ [days]"]
   model2=L"$\mathcal{H}_{PP}$"
   model3=L"$\mathcal{H}_{PPP}$"
   model=L"$\mathcal{H}_{PPPP}$"


  if obs=="fromEMB"
    name = string("OUTPUTS/EMBmc_table_",sigma,"s",nyear,"yrs.tex")
  else
    name = string("OUTPUTS/EVmc_table_",sigma,"s",nyear,"yrs.tex")
  end
  # open(name,"w") do io

    println("Model",'\t',model2,'\t',model3,'\t',model,"\\")
    println("BIC",'\t',BIC2,'\t',BIC3,'\t',BIC," \\")
 		println("χ^2",'\t',chi2,'\t',chi3,'\t',chi,"\\")
    println("reduced χ^2",'\t',reduced_chi2,'\t',reduced_chi3,'\t',reduced_chi,"\\")
  for i=1:length(avg)
    println(avg[i]," -",low[i]," +",high[i])
  end
  #   if i<=10
  #     println(io,parname[i],'\t',avg2[i],"_{-",avg2[i]-low2[i],"}^{+",high2[i]-avg2[i],"} & ",'\t',avg3[i],"_{-",avg3[i]-low3[i],"}^{+",high3[i]-avg3[i],"} & ",'\t',avg[i],"_{-",avg[i]-low[i],"}^{+",high[i]-avg[i],"} \\")
  #   end 
  #   if i >=11 && i <= 15
  #     println(io,parname[i],'\t','\t','\t',avg3[i],"_{-",avg3[i]-low3[i],"}^{+",high3[i]-avg3[i],"} & ",'\t',avg[i],"_{-",avg[i]-low[i],"}^{+",high[i]-avg[i],"} \\")
  #   end
  #   if i >=16 && i <= 20
  #     println(io,parname[i],'\t','\t','\t','\t','\t',avg[i],"_{-",avg[i]-low[i],"}^{+",high[i]-avg[i],"} \\")
  #   end
  # end
  # println(io,parname[end],'\t',avg2[end],"_{-",avg2[end]-low2[end],"}^{+",high2[end]-avg2[end],"} & ",'\t',avg3[end],"_{-",avg3[end]-low3[end],"}^{+",high3[end]-avg3[end],"} & ",'\t',avg[end],"_{-",avg[end]-low[end],"}^{+",high[end]-avg[end],"} \\")
  # end

  return 
end
