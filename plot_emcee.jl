using PyPlot

function plot_emcee(mcmc,pname,include_moon::Bool=false)
    #par_mcmc lprob_mcmc param nwalkers nsteps accept iburn indepsamples
  par_mcmc= mcmc["par_mcmc"]
  lprob_mcmc = mcmc["lprob_mcmc"]
  nwalkers = mcmc["nwalkers"]
  nsteps = mcmc["nsteps"]
  accept = mcmc["accept"]
  iburn = mcmc["iburn"]
  indepsamples = mcmc["indepsamples"]

  parname = [L"$μ_1$ [$M_{⋆}$]",L"$P_1$ [days]",L"$t_{0,1}$",L"$e_1 cos(ω_1)$",L"$e_1 sin(ω_1)$",
    L"$μ_2$ [$M_{⋆}$]",L"$P_2$ [days]",L"$t_{0,2}$",L"$e_2 cos(ω_2)$",L"$e_2 sin(ω_2)$",
    L"$μ_3$ [$M_{⋆}$]",L"$P_3$ [days]",L"$t_{0,3}$",L"$e_3 cos(ω_3)$",L"$e_3 sin(ω_3)$"]

  if string(pname) == "venus"
    figsize=(8,6)
    for i=1:5
    subplot(3,2,i)
    for j=1:nwalkers 
    plot(par_mcmc[j,iburn:nsteps,i])
    ylabel(parname[i])
    end
    end
    subplot(3,2,6)
    for j=1:nwalkers
      plot(lprob_mcmc[j,iburn:nsteps])  
      ylabel(L"$log_{10} Prob$")
    end
    tight_layout()
    name = string("IMAGES/MCMCstepsp1.png")
    # savefig(name)
    @show()
  elseif string(pname) == "earth"
    figsize=(8,6)
    for i=1:5
    subplot(3,2,i)
    for j=1:nwalkers 
    plot(par_mcmc[j,iburn:nsteps,i+5])
    ylabel(parname[i+5])
    end
    end
    subplot(3,2,6)
    for j=1:nwalkers
      plot(lprob_mcmc[j,iburn:nsteps])  
      ylabel(L"$log_{10} Prob$")
    end
    tight_layout()
    name = string("IMAGES/MCMCstepsp2.png")
    # savefig(name)
    @show()
  elseif string(pname) == "jupiter"
    figsize=(8,6)
    for i=1:5
    subplot(3,2,i)
    for j=1:nwalkers 
    plot(par_mcmc[j,iburn:nsteps,i+10])
    ylabel(parname[i+10])
    end
    end
    subplot(3,2,6)
    for j=1:nwalkers
      plot(lprob_mcmc[j,iburn:nsteps])  
      ylabel(L"$log_{10} Prob$")
    end
    tight_layout()
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
    plot(par_mcmc[j,iburn:nsteps,12])#,exp.(lprob_mcmc[j,1:nsteps]) .- maximum(lprob_mcmc[j,1:nsteps]))
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