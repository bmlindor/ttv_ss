using PyPlot

function plot_emcee(mcmc,pname,include_moon::Bool=false)
    #par_mcmc lprob_mcmc param nwalkers nsteps accept iburn indepsamples
  par_mcmc= mcmc["par_mcmc"]
  nwalkers = mcmc["nwalkers"]
  nsteps = mcmc["nsteps"]
  accept = mcmc["accept"]
  iburn = mcmc["iburn"]
  indepsamples = mcmc["indepsamples"]

  parname = ["mu_1","P_1","t01","e1 cos(om1)","e1 sin(om1)",
    "mu_2","P_2","t02","e2 cos(om2)","e2 sin(om2)",
    "mu_3","P_3","t03","e3 cos(om3)","e3 sin(om3)"]

  if string(pname) == "venus"
    figsize=(9,5)
    for i=1:5
      subplot(5,1,i)
      for j=1:nwalkers 
      plot(par_mcmc[j,1:nsteps,i])
      ylabel(parname[i])
      end
      tight_layout()
    end
    name = string("IMAGES/MCMCstepsp1.png")
    # savefig(name)
    @show()
  elseif string(pname) == "earth"
    figsize=(9,5)
    for i=1:5
      subplot(5,1,i)
      for j=1:nwalkers 
      plot(par_mcmc[j,1:nsteps,i+5])
      ylabel(parname[i+5])
      end
      tight_layout()
    end
    name = string("IMAGES/MCMCstepsp2.png")
    # savefig(name)
    @show()
  elseif string(pname) == "jup"
    figsize=(9,5)
    for i=1:5
      subplot(5,1,i)
      for j=1:nwalkers 
      plot(par_mcmc[j,1:nsteps,i+10])
      ylabel(parname[i+10])
      end
      tight_layout()
    end
    name = string("IMAGES/MCMCstepsp3.png")
    # savefig(name)
    @show()
  end
  if include_moon
    figsize=(5,3)
    for i=1:3
    subplot(3,1,i)
    for j=1:nwalkers 
    plot(par_mcmc[j,1:nsteps,i+15])
    ylabel(parname[i+15])
    end
    # tight_layout()
    end
    name = string("IMAGES/MCMCstepsmoon.png")
    # savefig(name)
  end
end

function plot_trace(mcmc)
  par_mcmc= mcmc["par_mcmc"]
  lprob_mcmc = mcmc["lprob_mcmc"]
  nwalkers = mcmc["nwalkers"]
  nsteps = mcmc["nsteps"]
  accept = mcmc["accept"]
  iburn = mcmc["iburn"]
  indepsamples = mcmc["indepsamples"]
# exp.((wide["lprob_p3"] .-maximum(p_33["lprob_p3"]))),12
  figsize=(8,6)
  subplot(2,1,1)
  for j=1:nwalkers
    plot(par_mcmc[j,iburn:nsteps,11])#,exp.(lprob_mcmc[j,1:nsteps]) .- maximum(lprob_mcmc[j,1:nsteps]))
    ylabel(L"$Period$")
  end
  subplot(2,1,2)
    for j=1:nwalkers
    plot(lprob_mcmc[j,iburn:nsteps])  
    ylabel(L"$logProb$")
    # xlim(0,10)
  end
end 
# Make plot of MCMC parameters after burn-in
function plot_parameters(xvalue,yvalue)
  par_mcmc= mcmc["par_mcmc"]
  nwalkers = mcmc["nwalkers"]
  nsteps = mcmc["nsteps"]
  accept = mcmc["accept"]
  iburn = mcmc["iburn"]
  indepsamples = mcmc["indepsamples"]



# figsize=(8,6)
# for i=2:nparam
#   for j=1:i-1
#     scatter(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]),vec(par_mcmc[1:nwalkers,iburn:nsteps,j]))
#     xlabel(pname[i])
#     ylabel(pname[j])
#   end
# end
# name = string("IMAGES/MCMCparams",label,".png")
end