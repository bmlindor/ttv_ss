include("bounds.jl")
include("CGS.jl")
include("misc.jl")
using TTVFaster,DelimitedFiles,JLD2
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
      if param[(iplanet-1)*5+1] < 0
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
  println("Independent Sample Size: ",indepsamples)

	# Calculate the auto-correlation function for each parameter:
	#acf=zeros(nparam)
	#for i=1:nparam
		#acf=autocorrelation(par_mcmc[:,iburn:nsteps,i],, )
	#end
	#println("Auto-Cor Function: ",acf)

  # Find mean and standard deviation of posteriors after burn-in:
	mean_posteriors=[mean(par_mcmc[:,iburn:nsteps,i]) for i=1:nparam]
	println("Mean posteriors after burn-in: ",mean_posteriors)
  sigtot=[sqrt((mean(vec(par_mcmc[:,iburn:nsteps,end])).*3600*24)^2 + (mean(sigtt).*3600*24)^2) ]
	println("Total timing uncertainty of: ",sigtot," secs.")


  #mean_mp = [mean(vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+1])).*CGS.MSUN/CGS.MEARTH for iplanet=1:nplanet]
  #mp_errs = [std(vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+1])).*abs(CGS.MSUN/CGS.MEARTH) for iplanet=1:nplanet]

  #mean_ecc=[mean(sqrt.(vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+4]).^2 .+ vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+5]).^2)) for iplanet=1:nplanet]
	#ecc_errs= [sqrt((vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+4]).^2 .* (std(vec((par_mcmc[:,iburn:nsteps,(iplanet-1)*5+4])).^2))  / (vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+4]).^2 .+ vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+5]).^2)) .+ (vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+5]).^2 .* (std(vec((par_mcmc[:,iburn:nsteps,(iplanet-1)*5+5])).^2))  / (vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+4]).^2 .+ vec(par_mcmc[:,iburn:nsteps,(iplanet-1)*5+5]).^2))) for iplanet=1:nplanet]


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

# function mc_BIC(tt0,nplanet,ntrans,par_mcmc,tt,sigtt,jmax,EM)
  # params=[median(par_mcmc[:,:,i]) for i=1:length(par_mcmc[1,end,:])]
 
# end
# Retrieve MCMC results after burn-in
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
  println("Added moon to model? ",include_moon)
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
  weight=ones(nt1+nt2)./ sigtt.^2 
  nparam=length(par_mcmc[1,end,:])
  # println("From fit: ",param)
  println("           Posterior Parameters from ",mcfile)
  # println(" lprob: ",lprob)
  avg=zeros(nparam)
  minus1sig=zeros(nparam)
  plus1sig=zeros(nparam)
  for i=1:nparam
   avg[i],minus1sig[i],plus1sig[i]=quantile(vec(par_mcmc[:,iburn:end,i]),[0.5,0.1587,0.8413])
   # println(pname[i]," = ",avg[i]," + ",abs(plus1sig[i]-avg[i])," _ ",abs(avg[i]-minus1sig[i]))
  end
  function chi_mcmc(tt0,nplanet,ntrans,par_mcmc,tt,sigtt,jmax,EM)
    chisq = 0.0  
    tt_model = TTVFaster.ttv_wrapper(tt0,nplanet,ntrans,par_mcmc[1:end-1],jmax,EM) 
    for j=1:length(tt)
      chisq += (tt[j]-tt_model[j])^2 / (sigtt[j]^2 + par_mcmc[end]^2)
    end
    return chisq
  end
  println("median posteriors: ",avg)
  mean_posteriors=[mean(par_mcmc[:,iburn:nsteps,i]) for i=1:nparam]
  println("mean posteriors: ",mean_posteriors)
  chi2_avg = chi_mcmc(tt0,nplanet,ntrans,avg,tt,sigtt,jmax,EM)
  chi2 = chi_mcmc(tt0,nplanet,ntrans,mean_posteriors,tt,sigtt,jmax,EM)
  N=length(tt0) ; k=nparam
  println(" χ^2 from median: ",chi2_avg," from mean: ",chi2)
  mc_BIC(chi2,k,N)=chi2 + k*log(N)
  println("BIC from median: ",mc_BIC(chi2_avg,k,N))
  println("BIC from mean: ",mc_BIC(chi2,k,N))

    # Percentage of walkers where diff. between median and quantile value is >100
    # bad_walk=0
    # for i in 1:nwalkers
    #    med,sig1,quant=quantile!(jldmc["par_mcmc"][i,jldmc["iburn"]:end,12],[0.5,0.68,0.9])
    #    if abs(quant-med)>100
    #    bad_walk+=1
    #    end
    # end
    # println("Bad walkers: ",bad_walk/nwalkers)

  masses=[median(vec(par_mcmc[:,iburn:end,i-4])) for i in 1:length(param) if i%5==0] .*CGS.MSUN/CGS.MEARTH
  mass_errs=[std(vec(par_mcmc[:,iburn:end,i-4])) for i in 1:length(param) if i%5==0] .*CGS.MSUN/CGS.MEARTH
  periods=[mean(vec(jldmc["par_mcmc"][:,iburn:end,i-3])) for i in 1:length(param) if i%5==0] 
  per_errs=[std(vec(jldmc["par_mcmc"][:,iburn:end,i-3])) for i in 1:length(param) if i%5==0] 
  ecc=[median(sqrt.(vec(par_mcmc[:,iburn:nsteps,(i-1)]).^2 .+ vec(par_mcmc[:,iburn:nsteps,(i)]).^2)) for i in 1:length(param) if i%5==0]

  #   ecc_err1 = calc_quad_errs(median(jldmc["par_mcmc"][:,iburn:end,4]),std(jldmc["par_mcmc"][:,iburn:end,4]),median(jldmc["par_mcmc"][:,iburn:end,5]),std(jldmc["par_mcmc"][:,iburn:end,5]))
  #   ecc_err2 = calc_quad_errs(median(jldmc["par_mcmc"][:,iburn:end,9]),std(jldmc["par_mcmc"][:,iburn:end,9]),median(jldmc["par_mcmc"][:,iburn:end,10]),std(jldmc["par_mcmc"][:,iburn:end,10]))
  #   ecc_errs = [ecc_err1,ecc_err2]

    # if grid_type_nplanet=="p3" || grid_type_nplanet=="p3moon"
  #     ecc_err3=calc_quad_errs(median(jldmc["par_mcmc"][:,iburn:end,14]),std(jldmc["par_mcmc"][:,iburn:end,14]),median(jldmc["par_mcmc"][:,iburn:end,15]),std(jldmc["par_mcmc"][:,iburn:end,15]))
  #     ecc_errs=[ecc_errs;ecc_err3]
    # elseif grid_type_nplanet=="p4" || grid_type_nplanet=="p3moonp4"
  #     ecc_err4=calc_quad_errs(median(jldmc["par_mcmc"][:,iburn:end,14]),std(jldmc["par_mcmc"][:,iburn:end,14]),median(jldmc["par_mcmc"][:,iburn:end,15]),std(jldmc["par_mcmc"][:,iburn:end,15]))
  #     ecc_err3=calc_quad_errs(median(jldmc["par_mcmc"][:,iburn:end,19]),std(jldmc["par_mcmc"][:,iburn:end,19]),median(jldmc["par_mcmc"][:,iburn:end,20]),std(jldmc["par_mcmc"][:,iburn:end,20]))
  #     ecc_errs=[ecc_errs;ecc_err4;ecc_err3]
    # end
  # if grid_type_nplanet=="p3moon" || grid_type_nplanet=="p3moonp4"
  #   tmax_errs=calc_quad_errs(median(jldmc["par_mcmc"][:,iburn:end,19]),std(jldmc["par_mcmc"][:,iburn:end,19]),median(jldmc["par_mcmc"][:,iburn:end,20]),std(jldmc["par_mcmc"][:,iburn:end,20]))
  # end

  sigsys=(median(vec(par_mcmc[:,iburn:end,end]))).* 3600*24
  sigsys_err=(std(vec(par_mcmc[:,iburn:end,end]))).* 3600*24
  sigtot=sqrt(sigsys^2 + sigma^2) 

  println("Retrieved values.")
  println("M_p[M⊕]=",masses," +/- ",mass_errs)
  println("Per [d]=",periods," +/- ",per_errs)
  # println("eccen. =",ecc," +/- ",ecc_errs)
  println("σsys[s]=",sigsys," +/- ",sigsys_err)
  println("σtot[s]=",sigtot)
    # return masses, mass_errs, periods, per_errs, sigtot
end