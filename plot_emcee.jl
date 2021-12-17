using PyPlot,Statistics,JLD2
# Plot MCMC mass vs period traces 
function plot_trace(sigma,nyear,sim,model,include_moon::Bool=false)
  if String(sim)=="EMB" && isfile(string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jl"))
    mcfile=string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  elseif isfile(string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  else 
    return  println("MCMC file for ",sim," with ",model," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  m=jldopen(String(mcfile),"r")
  par_mcmc= m["par_mcmc"]
  lprob_mcmc=m["lprob_mcmc"]
  nwalkers=m["nwalkers"]
  nsteps=m["nsteps"]
  accept=m["accept"]
  iburn=m["iburn"]
  indepsamples=m["indepsamples"]

  if include_moon || String(sim)=="p4"
    fig=figure(figsize=(8,6))
    subplot(421)
  else
    fig=figure(figsize=(7,5))
    subplot(321)
  end
  ax1=gca()
  for j=1:nwalkers
  ax1.plot(vec(par_mcmc[j,iburn:nsteps,1]).* CGS.MSUN/CGS.MEARTH)
  end
  ax1.set_ylabel(L"$μ_1$ [$M_{⋆}$]")
  if include_moon || String(sim)=="p4"
    subplot(422)
  else 
    subplot(322)
  end
  ax2=gca()
  for j=1:nwalkers
    ax2.plot(vec(par_mcmc[j,iburn:nsteps,2]))  
  end
  ax2.set_ylabel(L"$P_1$ [days]")
  if include_moon || String(sim)=="p4"
    subplot(423)
  else 
    subplot(323)
  end
  ax3=gca()
  for j=1:nwalkers
    ax3.plot(vec(par_mcmc[j,iburn:nsteps,6]).* CGS.MSUN/CGS.MEARTH)
  end
  ax3.set_ylabel(L"$μ_2$ [$M_{⋆}$]")
  if include_moon || String(sim)=="p4"
    subplot(424)
  else 
    subplot(324)
  end
  ax4=gca()
  for j=1:nwalkers
    ax4.plot(vec(par_mcmc[j,iburn:nsteps,7]))
  end
  ax4.set_ylabel(L"$P_2$ [days]")
  if include_moon || String(sim)=="p4"
    subplot(425)
  else
    subplot(325)
  end
    ax5=gca()
    for j=1:nwalkers
      ax5.plot(vec(par_mcmc[j,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH)
    end
    ax5.set_ylabel(L"$μ_3$ [$M_{⋆}$]")
  if include_moon || String(sim)=="p4"
    subplot(426)
  else
    subplot(326)
  end
    ax6=gca()
    for j=1:nwalkers
      ax6.plot(vec(par_mcmc[j,iburn:nsteps,12]))
    end
    ax6.set_ylabel(L"$P_3$ [days]")

  #   subplot(427)
  #   ax7=gca()
  #   for j=1:nwalkers
  #     ax7.plot(vec(par_mcmc[j,iburn:nsteps,16]).* CGS.MSUN/CGS.MEARTH)
  #   end
  #   ax7.set_ylabel(L"$μ_4$ [$M_{⋆}$]")
  #   subplot(428)
  #   ax8=gca()
  #   for j=1:nwalkers
  #     ax8.plot(vec(par_mcmc[j,iburn:nsteps,16]))
  #   end
  #   ax8.set_ylabel(L"$P_4$ [days]")
  # end
  # if include_moon 
  #   subplot(425)
  # else 
  #   subplot(325)
  # end
  # ax5=gca()
  #   for j=1:nwalkers
  #     ax5.plot(vec(par_mcmc[j,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH)
  #   end
  #   ax5.set_ylabel(L"$μ_3$ [$M_{⋆}$]")
  # if include_moon 
  #   subplot(426)
  # else 
  #   subplot(326)
  # end
  # ax6=gca()
  # for j=1:nwalkers
  #   ax6.plot(vec(par_mcmc[j,iburn:nsteps,12]))
  # end
  # ax6.set_ylabel(L"$P_3$ [days]")
  # if include_moon
  #   ax7=subplot(4,2,7)
  #   for j=1:nwalkers
  #     ax7.plot(vec(par_mcmc[j,iburn:nsteps,18]))
  #   end
  #   ax7.set_ylabel(L"$Δϕ$ [rad]")
  #   ax8=subplot(4,2,8)
  #   for j=1:nwalkers
  #     ax8.plot(sqrt.(vec(par_mcmc[j,iburn:nsteps,16]).^2 + vec(par_mcmc[j,iburn:nsteps,17]).^2))
  #   end
  #   ax8.set_ylabel(L"$t_{max}$ [days]")
  # end
  # # ax3.set_ylabel(L"$σ_{sys}^2$ [days]")
  # tight_layout()
  # title=string("IMAGES/traces/",sim,model,"traces-",sigma,"secs",nyear,"yrs.png")
  # savefig(title)
end 
# Plot MCMC traces of individual walkers
function plot_emcee(sigma,nyear,sim,model,include_moon::Bool=false)
  if String(sim)=="EMB" && isfile(string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  elseif isfile(string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  else 
    return  println("MCMC file for ",sim," with ",model," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  m=jldopen(String(mcfile),"r")
  par_mcmc= m["par_mcmc"]
  lprob_mcmc=m["lprob_mcmc"]
  nwalkers=m["nwalkers"]
  nsteps=m["nsteps"]
  accept=m["accept"]
  iburn=m["iburn"]
  indepsamples=m["indepsamples"]
  parname=[L"$μ_1$ [$M_{⋆}$]",L"$P_1$ [days]",L"$t_{0,1}$",L"$e_1 cos(ω_1)$",L"$e_1 sin(ω_1)$",
    L"$μ_2$ [$M_{⋆}$]",L"$P_2$ [days]",L"$t_{0,2}$",L"$e_2 cos(ω_2)$",L"$e_2 sin(ω_2)$",
    L"$μ_3$ [$M_{⋆}$]",L"$P_3$ [days]",L"$t_{0,3}$",L"$e_3 cos(ω_3)$",L"$e_3 sin(ω_3)$",
    L"$μ_4$ [$M_{⋆}$]",L"$P_4$ [days]",L"$t_{0,4}$",L"$e_4 cos(ω_4)$",L"$e_4 sin(ω_4)$",
    L"$t_{max} sin(ϕ_0)$",L"$t_{max} cos(ϕ_0)$",L"$Δϕ$ [rad]",L"$σ_{sys}^2$ [days]"]
  figure(figsize=(8,6))
  for i=1:5
  ax1=subplot(3,2,i)
  for j=1:nwalkers 
  ax1.plot(par_mcmc[j,iburn:nsteps,i])
  end
  ax1.set_ylabel(parname[i])
  end
  tight_layout()
  title=string("IMAGES/traces/",sim,model,"testVenus-",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  figure(figsize=(8,6))
  for i=1:5
  ax2=subplot(3,2,i)
  for j=1:nwalkers 
  ax2.plot(par_mcmc[j,iburn:nsteps,i+5])
  end
  ax2.set_ylabel(parname[i+5])
  end
  tight_layout()
  title=string("IMAGES/traces/",sim,model,"testEarth-",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  if String(model)=="p4"
    figure(figsize=(8,6))
    for i=1:5
    ax3=subplot(3,2,i)
    for j=1:nwalkers 
    ax3.plot(par_mcmc[j,iburn:nsteps,i+10])
    end
    ax3.set_ylabel(parname[i+10])
    end
    tight_layout()
    title=string("IMAGES/traces/",sim,model,"testMars-",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
    figure(figsize=(8,6))
    for i=1:5
    ax3=subplot(3,2,i)
    for j=1:nwalkers 
    ax3.plot(par_mcmc[j,iburn:nsteps,i+15])
    end
    ax3.set_ylabel(parname[i+15])
    end
    ax4=subplot(3,2,6)
    for j=1:nwalkers
      ax4.plot(par_mcmc[j,iburn:nsteps,end])
      ax4.set_ylabel(parname[end])
    end
    tight_layout()
    title=string("IMAGES/traces/",sim,model,"testJupiter-",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
  else
    figure(figsize=(8,6))
    for i=1:5
    ax3=subplot(3,2,i)
    for j=1:nwalkers 
    ax3.plot(par_mcmc[j,iburn:nsteps,i+10])
    end
    ax3.set_ylabel(parname[i+10])
    end
    ax4=subplot(3,2,6)
    for j=1:nwalkers
      ax4.plot(par_mcmc[j,iburn:nsteps,end])
      ax4.set_ylabel(parname[end])
    end
    tight_layout()
    title=string("IMAGES/traces/",sim,model,"testJupiter-",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
  end
  if include_moon
    figure(figsize=(7,5))
    for i=1:3
    ax5=subplot(2,2,i)
    for j=1:nwalkers 
    ax5.plot(par_mcmc[j,1:nsteps,i+15])
    end
    ax5.set_ylabel(parname[i+20])
    end
    # subplot(2,2,4)
    # for j=1:nwalkers
    # plot(lprob_mcmc[j,iburn:nsteps])  
    # ylabel(L"$logProb$")
    # end
    tight_layout()
    title=string("IMAGES/traces/",sim,model,"Moon-",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
  end
end

function plot_effects(sigma,nyear,sim,model)
  if String(sim)=="EMB" && isfile(string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  elseif isfile(string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  else 
    return  println("MCMC file for ",sim," with ",model," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  m=jldopen(String(mcfile),"r")
  par_mcmc=m["par_mcmc"]
  lprob_mcmc=m["lprob_mcmc"]
  nwalkers=m["nwalkers"]
  nsteps=m["nsteps"]
  accept=m["accept"]
  iburn=m["iburn"]
  plot()
end
# # figsize=(8,6)
# # for i=2:nparam
# #   for j=1:i-1
# #     scatter(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]),vec(par_mcmc[1:nwalkers,iburn:nsteps,j]))
# #     xlabel(pname[i])
# #     ylabel(pname[j])
# #   end
# # end
# # name=string("IMAGES/MCMCparams",label,".png")