  

# Run a Markov chain:
function MCMC() 
    errors = [1e-7,1e-5,1e-5,1e-2,1e-2,1e-7,1e-5,1e-5,1e-2,1e-2,1e-6,1e-1,1e-1,1e-2,1e-2]
    pname = ["mu_1","P_1","t01","e1 cos(om1)","e1 sin(om1)","mu_2","P_2","t02","e2 cos(om2)","e2 sin(om2)","mu_3","P_3","t03","e3 cos(om3)","e3 sin(om3)"]
    nwalkers = nparam * 3
    nsteps = 10000
    #nsteps = 100
    # Set up arrays to hold the results:
    par_mcmc = zeros(nwalkers,nsteps,nparam)
    chi_mcmc = zeros(nwalkers,nsteps)
    # Initialize walkers:
    par_trial = fit.param
    for j=1:nwalkers
    # Select from within uncertainties:
      chi_trial = 1e100
    # Only initiate models with reasonable chi-square values:
      while chi_trial > chi_best + 1000
        par_trial = fit.param + errors.*randn(nparam)
        # model = ttv_wrapper3(tt0,par_trial)
        model = ttv_wrapper(tt0, nplanet, ntrans, par_trial)
        chi_trial = sum(((tt-model)./sigtt).^2)
        println("chi_trial: ",chi_trial)
      end
      chi_mcmc[j,1]=chi_trial
      par_mcmc[j,1,:]=par_trial
      println("Success: ",par_trial,chi_trial)
    end
    # Initialize scale length & acceptance counter:
    ascale = 2.0
    accept = 0
    # Next, loop over steps in markov chain:
    for i=2:nsteps
      for j=1:nwalkers
        ipartner = j
    # Choose another walker to 'breed' a trial step with:
        while ipartner == j
          ipartner = ceil(Int,rand()*nwalkers)
        end
    # Now, choose a new set of parameters for the trial step:
        z=(rand()*(sqrt(ascale)-1.0/sqrt(ascale))+1.0/sqrt(ascale))^2
        par_trial=vec(z*par_mcmc[j,i-1,:]+(1.0-z)*par_mcmc[ipartner,i-1,:])
    # Compute model & chi-square:  
        # model_trial =ttv_wrapper3(tt0,par_trial)
        model_trial = ttv_wrapper(tt0, nplanet, ntrans, par_trial)
        chi_trial=sum(((tt-model_trial)./sigtt).^2)
    # Next, determine whether to accept this trial step:
        alp = z^(nparam-1)*exp(-0.5*(chi_trial - chi_mcmc[j,i-1]))
        if rand() < 0.0001
          println("Step: ",i," Walker: ",j," Chi-square: ",chi_trial," Prob: ",alp," Frac: ",accept/(mod(i-1,1000)*nwalkers+j))
        end
        if alp >= rand()
    # If step is accepted, add it to the chains!
          par_mcmc[j,i,:] = par_trial
          chi_mcmc[j,i] = chi_trial
          accept = accept + 1
        else
    # If step is rejected, then copy last step:
          par_mcmc[j,i,:] = par_mcmc[j,i-1,:]
          chi_mcmc[j,i] = chi_mcmc[j,i-1]
        end
      end
      if mod(i,1000) == 0
        println("Number of steps: ",i," acceptance rate: ",accept/(1000*nwalkers))
        accept = 0
      end
    end

    for i=1:nparam
      for j=1:nwalkers
        plot(vec(par_mcmc[j,1:nsteps,i]))
      end
      xlabel("MCMC step")
      ylabel(pname[i])
      println("Hit return to continue")
      read(stdin,Char)
      clf()
    end

    # Now, determine time of burn-in by calculating first time median is crossed:
    iburn = 0
    for i=1:nparam
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

    println("Burn-in ends: ",iburn)

    for i=2:nparam
      for j=1:i-1
        scatter(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]),vec(par_mcmc[1:nwalkers,iburn:nsteps,j]))
        xlabel(pname[i])
        ylabel(pname[j])
        println("Hit return to continue")
        read(stdin,Char)
        clf()
      end
    end
    return par_mcmc,chi_mcmc