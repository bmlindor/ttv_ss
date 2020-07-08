include("ttv_wrapper2.jl")
include("ttv_wrapper3.jl")
include("ttv_wrapper_fixp3.jl")
include("regress.jl")
include("chisquare3.jl")
include("chisquare2.jl")
include("compute_ttv.jl")
using LsqFit
using PyPlot
using Optim
using DelimitedFiles

function fit_mysteryplanet3()
    #=
     To do:
     1). Carry out a linear fit to the transit times.
     2). Write a wrapper to call ttv_nplanet.jl which then
         computes the chi-square of the fit.
     3). Do a fit to the inner two planets.
     4). Initialize a grid of periods & phases of the outer planet,
         & compute the best-fit at each with an optimizer.
     5). Write a markov chain to compute the best-fit parameters
         for the 3 planets.
    =#

    data1 = readdlm("ttv_venus.txt")
    tt1 = vec(data1[:,1])
    nt1 = length(tt1)
    data2 = readdlm("ttv_earth.txt")
    tt2 = vec(data2[:,1])
    nt2 = length(tt2)

    # Okay, let's do a linear fit to the transit times (first column):
      #linear fit that we already did
    x = zeros(2,nt1)
    x[1,1:nt1] .= 1.0
    x[2,1:nt1]=range(0,stop=nt1-1,length=nt1)
    sigtt1 = ones(nt1).* 30 ./ 24 ./3600
    coeff1, covcoeff1 = regress(x,tt1,sigtt1)
    t01 = coeff1[1]; per1 = coeff1[2]

    x = zeros(2,nt2)
    x[1,1:nt2].=1.0
    x[2,1:nt2]=range(0,stop=nt2-1,length=nt2)
    sigtt2 = ones(nt2).* 30 ./ 24 ./3600
    coeff2, covcoeff2 = regress(x,tt2,sigtt2)
    sigtt=[sigtt1;sigtt2]
    t02 = coeff2[1]; per2 = coeff2[2]
    t1  = collect(t01 .+ per1 .* range(0,stop=nt1-1,length=nt1)) #best fit linear transit times w/o ttvs
    t2  = collect(t02 .+ per2 .* range(0,stop=nt2-1,length=nt2))
    # Best-fit linear transit times:
    tt0 = [t1;t2]
    weight = ones(nt1+nt2) #assigns each data point stat weight d.t. noise
    # Actual transit times:
    tt=[tt1;tt2]

    # Okay, now let's do a 2-planet fit:
        #mass ratio, period, initial transit time, e*cos(omega), e*sin(omega)
    param = [3e-6,per1,t01,0.01,0.01,3e-6,per2,t02,0.01,0.01] #initial params defined for both planets
    println("Initial parameters: ",param)
    #model = ttv_wrapper2(tt0,param)

        #sets up data structure to hold planet properties, passed to TTVFaster
    jmax = 5
    data=param
    n1 = length(tt1);
    # global n1
    n2 = length(tt2);
    # global n2
    p1=TTVFaster.Planet_plane_hk(data[1],data[2],data[3],data[4],data[ 5])
    p2=TTVFaster.Planet_plane_hk(data[6],data[7],data[8],data[9],data[10])
    time1 = collect(p1.trans0 .+ range(0,stop=n1-1,length=n1) .* p1.period)
    time2 = collect(p2.trans0 .+ range(0,stop=n2-1,length=n2) .* p2.period)
    # Initialize the computation of the Laplace coefficients:
    ttv1 = zeros(n1)
    ttv2 = zeros(n2)

    dummy=TTVFaster.compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2) #first call to TTVFaster w/o optim

    scatter(time1,tt1.-t1)
    plot(time1,ttv1)
    scatter(time2,tt2.-t2,color="green")
    plot(time2,ttv2)

        #globals are used by TTVFaster, actual observed data held fixed in optim
    global tt0
    global tt
    global sigtt
    p3_cur = 11.86*365.25 #jupiter period in days, initial value
    global p3_cur

    #res = optimize(chisquare2, param, method = :l_bfgs, iterations = 21)
    res = optimize(chisquare2, param) #optimizes 2 planet fit

    #fit2 = curve_fit(ttv_wrapper2,tt0,tt,weight,param; show_trace=true)

    #println("Finished 2-planet fit: ",param)

    # Now, let's add the 3rd planet:
    #p3 = 11.86*365.25
    # Number of periods of the 3rd planet to search over:
    p3in  = 4000
    p3out = 4600
    np3 = 10
    #p3in  = 500 #larger range in period
    #p3out = 10000
    #np3 = 1000
    #p3in  = 11*365
    #p3out = 13*365
    #p3 = linspace(p3in,p3out,np3)
    p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
    nphase = 10
    chi_p3 = zeros(np3)
    nparam = 15
    param_p3 = zeros(nparam,np3)
    chi_best = 1e100 #global best fit
    pbest = zeros(nparam)
    for j=1:np3
      phase = p3[j]*range(0,stop=1,length=nphase) #searches over period of jupiter
      chi_phase = zeros(nphase)
      chi_p3[j] = 1e100
      for i=1:nphase #loops over jupiter phases
        param_tmp = [1e-3,phase[i],0.01,0.01] 
        param3 = [param;param_tmp] #concatenate 2 planet model to 3 planet model params
        p3_cur = p3[j] #sets jupiter period to global value
        fit = curve_fit(ttv_wrapper_fixp3,tt0,tt,weight,param3) #optimizes fit w/ 3 planet model
        ttmodel=ttv_wrapper_fixp3(tt0,fit.param)
        chi_phase[i]= sum((tt-ttmodel).^2 ./sigtt.^2)
        if chi_phase[i] < chi_best # check that best fit for period is better than global best fit
          chi_best = chi_phase[i]
          pbest = [fit.param[1:11];p3_cur;fit.param[12:14]]
        end
        if chi_phase[i] < chi_p3[j] # checks best fit over all phases of jupiter for this particular period
          chi_p3[j] = chi_phase[i]
          param_p3[1:nparam,j] =  [fit.param[1:11];p3_cur;fit.param[12:14]]
        end
      end
      println("Period: ",p3[j]," chi: ",chi_p3[j]," Param: ",vec(param_p3[1:nparam,j]))
    end
    clf()
    plot(p3/365.25,exp(-0.5*(chi_p3-minimum(chi_p3)))) #to show that max likelihood peaks at actual period
    xlabel("Period of planet 3 [years]")
    ylabel("Likelihood")
    println("Hit return to continue")
    read(STDIN,Char)
    clf()

    #
    #println("Calling ttv_wrapper3")
    #ttmodel=ttv_wrapper3(tt0,param3)
    #println("Finished ttv_wrapper3")
    #res = optimize(chisquare3, param3, method = :l_bfgs, iterations = 21)
    #  res = optimize(chisquare3, param3, method = :l_bfgs)
    #  ttmodel=ttv_wrapper3(tt0,param3)
    println("Best-fit parameters: ",pbest)
    fit = curve_fit(ttv_wrapper3,tt0,tt,weight,pbest)
    ttmodel=ttv_wrapper3(tt0,pbest)
    chi_best= sum((tt-ttmodel).^2 ./sigtt.^2)
    println("Minimum: ",chi_best," Param: ",fit.param)
    println(fit.param)

    scatter(time1,tt1.-t1)
    plot(time1,ttmodel[1:n1].-t1)
    scatter(time2,tt2.-t2,color="green")
    plot(time2,ttmodel[n1+1:n1+n2].-t2)

    println("Hit return to continue")
    read(STDIN,Char)
    clf()

    #println(fit2.param)
    #end

    # Run a Markov chain:
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
        model = ttv_wrapper3(tt0,par_trial)
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
        model_trial =ttv_wrapper3(tt0,par_trial)
        chi_trial=sum(((tt-model_trial)./sigtt).^2)
    # Next, determine whether to accept this trial step:
        alp = z^(nparam-1)*exp(-0.5*(chi_trial - chi_mcmc[j,i-1]))
        if rand() < 0.0001
          println("Step: ",i," Walker: ",j," Chi-square: ",chi_trial," Prob: ",alp," Frac: ",accept/(mod(i-1,1000)*nwalkers+j))
        end
        if alp >= rand()
    # If step is accepted, add it to the chains!
          par_mcmc[j,i,:] = par_trial
          chi_mcmc[j,i,:] = chi_trial
          accept = accept + 1
        else
    # If step is rejected, then copy last step:
          par_mcmc[j,i,:] = par_mcmc[j,i-1,:]
          chi_mcmc[j,i,:] = chi_mcmc[j,i-1]
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
      read(STDIN,Char)
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
        read(STDIN,Char)
        clf()
      end
    end
    return par_mcmc,chi_mcmc
end
