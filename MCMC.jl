if !@isdefined(TTVFaster)
    include("TTVFaster/src/TTVFaster.jl")
    using Main.TTVFaster
end
import Main.TTVFaster.ttv_wrapper
import Main.TTVFaster.chisquare
include("bounds.jl")
using DelimitedFiles,JLD2,Statistics,MCMCDiagnostics

# Run a Markov chain:
function MCMC(foutput::String,param::Array{Float64,1},lprob_best::Float64,
  nsteps::Int64,nwalkers::Int64,nplanet::Int64,ntrans::Array{Int64,1},
  tt0::Array{Float64,1},tt::Array{Float64,1},sigtt::Array{Float64,1},
  use_sigsys::Bool,EMB::Bool) 

  nparam = length(param)
  jmax = 5
  errors = [1e-7,1e-5,1e-5,1e-2,1e-2,
            1e-7,1e-5,1e-5,1e-2,1e-2,
            1e-6,1e-1,1e-1,1e-2,1e-2]
  pname = ["mu_1","P_1","t01","e1 cos(om1)","e1 sin(om1)",
          "mu_2","P_2","t02","e2 cos(om2)","e2 sin(om2)",
          "mu_3","P_3","t03","e3 cos(om3)","e3 sin(om3)"]
  if EMB
    errors = errors
    pname = pname
  else
    moon_errors = [1e-4,1e-4,1e-5]
    moon_name = ["tmax sin(phi0)","tmax cos(phi0)","deltaphi"]
    append!(errors,moon_errors)
    append!(pname,moon_name)
  end
  # Initialize walkers:
  if use_sigsys
    par_trial = [param[1:end];1e-8]
    nsize = nparam+1
  else
    par_trial = param[1:end]
    nsize = nparam
  end
  @assert (nwalkers >= 40)  # nwalkers = nparam * 3
  # Set up arrays to hold the results:
  par_mcmc = zeros(nwalkers,nsteps,nsize)
  lprob_mcmc = zeros(nwalkers,nsteps)

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
    if !EMB 
      # Plase priors on deltaphi and account for aliasing:
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
      # sigsys priors:
      sigsys = param[end]
      if sigsys < 0 
        lprior += -Inf
      end
    end
    return lprior
  end

  Nobs = sum(ntrans) - 2
  for j=1:nwalkers
  # Select from within uncertainties:
    lprob_trial = -1e100
  # Only initiate models with reasonable chi-square values:
    while lprob_trial < lprob_best - 1000 #since logProb, maybe 1000 is too large
      par_trial[1:nparam] .= param .+ errors .* randn(nparam) 
      if use_sigsys
        par_trial[nparam+1] = 1e-8 .* abs(randn())
      end
      lprob_trial = calc_lprior(par_trial)
      # println("Calculated Log Prior: ",lprob_trial)
      if lprob_trial > -Inf
        # model = ttv_wrapper3(tt0,par_trial)# 
        # if EMB
        #   model = ttv_wrapper(tt0,nplanet,ntrans,par_trial[1:nparam],jmax,true)
        # else 
          model = ttv_wrapper(tt0,nplanet,ntrans,par_trial[1:nparam],jmax,EMB)
        # end
        # println(length(tt)," ",length(model))
        if use_sigsys
          # ll = (-0.5 * sum((tt-model).^2 ./(sigtt.^2 .+ sigsys.^2) .+ log.(sigtt.^2 .+ sigsys.^2)))
         lprob_trial += (-0.5 * sum((tt-model).^2 ./(sigtt.^2 .+ par_trial[end]) .+ log.(sigtt.^2 .+ par_trial[end])))
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
  # Compute model & chi-square:  
      # model_trial =ttv_wrapper3(tt0,par_trial)
      lprob_trial = calc_lprior(par_trial)  
      if lprob_trial > -Inf
        # if EMB
        #   model_trial = ttv_wrapper(tt0,nplanet,ntrans,par_trial[1:nparam],jmax,true)
        # else 
          model_trial = ttv_wrapper(tt0,nplanet,ntrans,par_trial[1:nparam],jmax,EMB)
        # end
        # model_trial = ttv_wrapper(tt0,nplanet,ntrans,par_trial[1:nparam],false)
        # ll = -.5 *sum((tt-model_trial).^2 ./(sigtt.^2 .+ par_trial[end]) .+ log.(sigtt.^2 .+ par_trial[end]))
        # ll =  log(sum((tt-model_trial).^2 ./sigtt.^2))
        # lprob_trial = calc_lprior(par_trial) + ll*(1 - Nobs/2) 
        if use_sigsys
         lprob_trial += (-0.5 * sum((tt-model_trial).^2 ./(sigtt.^2 .+ par_trial[end]) .+ log.(sigtt.^2 .+ par_trial[end])))
        else
          lprob_trial += log(sum((tt-model_trial).^2 ./sigtt.^2))*(1 - Nobs/2) 
        end 
      end   
  # Next,determine whether to accept this trial step:
      alp = z^(nsize-1)*exp((lprob_trial - lprob_mcmc[j,i-1])) #prob ratio between current trial and previous prob
      if rand() < 0.0001
        println("Step: ",i," Walker: ",j," Trial Log Prob: " ,lprob_trial," Prob: ",alp," Frac: ",accept/(mod(i-1,1000)*nwalkers+j))
      end
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
    if mod(i,1000) == 0
      println("Number of steps: ",i," Acceptance Rate: ",accept/(1000*nwalkers))
      accept = 0
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

  # Now,calculate the minimum number of effective samples over all parameters:
  samplesize = zeros(nparam)
  for j=1:nwalkers
    for i=1:nparam
      samplesize[i] += effective_sample_size(par_mcmc[j,:,i])
    end
  end
  indepsamples = minimum(samplesize)
  println("Independent Sample Size: ",indepsamples)
  
  mcmcfile = string(foutput)
  @save mcmcfile par_mcmc lprob_mcmc param nwalkers nsteps accept iburn indepsamples
  return lprob_mcmc,par_mcmc #, param, nwalkers, nsteps, accept, iburn, indepsamples
end