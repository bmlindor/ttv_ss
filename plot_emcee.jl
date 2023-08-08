using PyPlot,Statistics,JLD2
include("CGS.jl")
# Plot MCMC mass vs period traces 
function plot_trace(sigma::Real,nyear::Real,grid_type_nplanet::String,case_num,include_moon::Bool=false)
  if case_num==1 #&& isfile(string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jl"))
    mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    case_label="Case 1"
    title="Posteriors from Venus + EMB TTVs"
  elseif case_num==2 #&& isfile(string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    case_label="Case 2"
    title="Posteriors from Venus + Earth TTVs"
  elseif case_num==3
    mcfile=string("MCMC/no_noise/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    title="Search from Venus + Earth TTVs -- no noise"
  else
    return  println("MCMC file for ",grid_type_nplanet," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  m=jldopen(String(mcfile),"r")
  par_mcmc= m["par_mcmc"]
  lprob_mcmc=m["lprob_mcmc"]
  nwalkers=m["nwalkers"]
  nsteps=m["nsteps"]
  pname=m["pname"]
  # accept=m["accept"]
  iburn=m["iburn"]
  indepsamples=m["indepsamples"]
  sigsys=(mean(vec(m["par_mcmc"][:,iburn:end,end]))).* 3600*24
  sigsys_err=(std(vec(m["par_mcmc"][:,iburn:end,end]))).* 3600*24
  sigtot=sqrt(sigsys^2 + sigma^2) 
  fig,ax=subplots(figsize=(8,6))
  suptitle(string(title,'\n',nyear," yr span ",L"$\sigma_{obs}=$",sigma,"s ",L"$\sigma_{sys}=$",round(sigsys,sigdigits=4),"s"))
  subplots_adjust(hspace=0.01,wspace=0.15)
  # text(0.5,0.5,string()
  function main_traces()
    # ax=fig.add_axes([0.5,0.2,0.3,0.5])
    ax1=subplot(421) #trace plot
    for j=1:nwalkers
      ax1.plot(vec(par_mcmc[j,iburn:nsteps,1]).* CGS.MSUN/CGS.MEARTH)
    end
    ax1.tick_params(labelbottom=false)
    ax1.set_ylabel(L"$μ_1$ [$M_{⋆}$]")
    ax2=subplot(422)
    for j=1:nwalkers
      ax2.plot(vec(par_mcmc[j,iburn:nsteps,2]))  
    end
    ax2.tick_params(labelbottom=false,labelleft=false,labelright=true,right=true)
    ax2.set_ylabel(L"$P_1$ [days]")
    ax3=subplot(423)
    for j=1:nwalkers
      ax3.plot(vec(par_mcmc[j,iburn:nsteps,6]).* CGS.MSUN/CGS.MEARTH)
    end
    ax3.tick_params(labelbottom=false)
    ax3.set_ylabel(L"$μ_2$ [$M_{⋆}$]")
    ax4=subplot(424)
    for j=1:nwalkers
      ax4.plot(vec(par_mcmc[j,iburn:nsteps,7]))
    end
    ax4.tick_params(labelbottom=false,labelleft=false,labelright=true,right=true)
    ax4.set_ylabel(L"$P_2$ [days]")
    if length(pname)>19
      ax5=subplot(425)
      for j=1:nwalkers
        ax5.plot(vec(par_mcmc[j,iburn:nsteps,16]).* CGS.MSUN/CGS.MEARTH)
      end
      ax5.set_ylabel(L"$μ_3$ [$M_{⋆}$]")
      ax6=subplot(426)
      for j=1:nwalkers
        ax6.plot(vec(par_mcmc[j,iburn:nsteps,17]))
      end
      ax6.set_ylabel(L"$P_3$ [days]")
      ax5.tick_params(labelbottom=false)
      ax6.tick_params(labelbottom=false,labelleft=false,labelright=true,right=true)
      ax7=subplot(427)
      for j=1:nwalkers
        ax7.plot(vec(par_mcmc[j,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH)
      end
      ax7.set_ylabel(L"$μ_4$ [$M_{⋆}$]")
      ax7.set_xlabel("No. Steps")
      ax8=subplot(428)
      for j=1:nwalkers
        plot(vec(par_mcmc[j,iburn:nsteps,12]))
      end
      ax8.set_xlabel("No. Steps")
      ax8.set_ylabel(L"$P_4$ [days]")
      ax8.tick_params(labelleft=false,labelright=true,right=true)
    else
      ax5=subplot(425)
      for j=1:nwalkers
        ax5.plot(vec(par_mcmc[j,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH)
      end
      ax5.set_ylabel(L"$μ_3$ [$M_{⋆}$]")
      ax6=subplot(426)
      for j=1:nwalkers
        ax6.plot(vec(par_mcmc[j,iburn:nsteps,12]))
      end
      ax6.set_ylabel(L"$P_3$ [days]")
      ax5.set_xlabel("No. Steps")
      ax6.set_xlabel("No. Steps")
      ax6.tick_params(labelbottom=true,labelleft=false,labelright=true,right=true)
    end
  end
##plotting vs prob
  function jup_traces()
    ax1=subplot(221) # planet 4 probability traces
    suptitle(string(title,'\n'," [",nyear," yr span] ",L"$\sigma_{obs}=$",sigma," sec ",L"$\sigma_{sys}=$",round(sigsys,sigdigits=3)," sec"))
      for j=1:nwalkers
      plot(vec(par_mcmc[j,iburn:nsteps,end]).* 3600*24,lprob_mcmc[j,iburn:nsteps])
      xlabel(L"$\sigma_{sys}$ [sec]")
        # plot(vec(lprob_mcmc[j,iburn:nsteps]),label=string("walker=",j))
      end
      subplot(222)
      for j=1:nwalkers
      plot(vec(par_mcmc[j,iburn:nsteps,16]).* CGS.MSUN/CGS.MEARTH,lprob_mcmc[j,iburn:nsteps])
      xlabel(L"Mass [M$\oplus$]")
      end
      subplot(223)
      for j=1:nwalkers
      plot(vec(par_mcmc[j,iburn:nsteps,17]),lprob_mcmc[j,iburn:nsteps])
      xlabel("Period [days]")
      end
      subplot(224)
      for j=1:nwalkers
      ecc=vec(sqrt.(par_mcmc[j,iburn:nsteps,19].^2 .+ par_mcmc[j,iburn:nsteps,20].^2))
      plot(ecc,lprob_mcmc[j,iburn:nsteps])
      xlabel("Eccentricity")
      end
    tight_layout()
    name=string("IMAGES/trace/case",case_num,"Jup",grid_type_nplanet,sigma,"secs",nyear,"yrs.png")
  end
##plotting vs prob
  function mars_traces()
    ax1=subplot(221) #planet 3 probabbility traces
    suptitle(string(title,'\n'," [",nyear," yr span] ",L"$\sigma_{obs}=$",sigma," sec ",L"$\sigma_{sys}=$",round(sigsys,sigdigits=3)," sec"))
      for j=1:nwalkers
      plot(vec(par_mcmc[j,iburn:nsteps,end]).* 3600*24,lprob_mcmc[j,iburn:nsteps])
      xlabel(L"$\sigma_{sys}$ [set]")
      end
      subplot(222)
      for j=1:nwalkers
      plot(vec(par_mcmc[j,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH,lprob_mcmc[j,iburn:nsteps])
      xlabel(L"Mass [M$\oplus$]")
      end
      subplot(223)
      for j=1:nwalkers
        # if j=/= 56
      plot(vec(par_mcmc[j,iburn:nsteps,12]),lprob_mcmc[j,iburn:nsteps])
      end
      xlabel("Period [days]")
      subplot(224)
      for j=1:nwalkers
      ecc=vec(sqrt.(par_mcmc[j,iburn:nsteps,14].^2 .+ par_mcmc[j,iburn:nsteps,15].^2))
      plot(ecc,lprob_mcmc[j,iburn:nsteps])
      xlabel("Eccentricity")
      end
    name=string("IMAGES/trace/case",case_num,"Mar",grid_type_nplanet,sigma,"secs",nyear,"yrs.png")
  end
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
  # mars_traces()
  name=string("IMAGES/trace/case",case_num,"main",grid_type_nplanet,sigma,"s",nyear,"yrs.png")
  main_traces()
  # tight_layout()
  savefig(name)
    # clf()
end 
# Plot MCMC traces of individual walkers
function plot_emcee(sigma::Real,nyear::Real,sim::String,model::String,include_moon::Bool=false)
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
    L"$μ_5$ [$M_{⋆}$]",L"$P_5$ [days]",L"$t_{0,5}$",L"$e_4 cos(ω_5)$",L"$e_5 sin(ω_5)$",
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
  title=string("IMAGES/traces/",sim,model,"Venus-",sigma,"secs",nyear,"yrs.png")
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
  title=string("IMAGES/traces/",sim,model,"Earth-",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  if String(model)=="p5"
    figure(figsize=(8,6))
    for i=1:5
    ax3=subplot(3,2,i)
    for j=1:nwalkers 
    ax3.plot(par_mcmc[j,iburn:nsteps,i+20])
    end
    ax3.set_ylabel(parname[i+20])
    end
    tight_layout()
    title=string("IMAGES/traces/",sim,model,"Saturn-",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
    figure(figsize=(8,6))
    for i=1:5
    ax3=subplot(3,2,i)
    for j=1:nwalkers 
    ax3.plot(par_mcmc[j,iburn:nsteps,i+10])
    end
    ax3.set_ylabel(parname[i+10])
    end
    tight_layout()
    title=string("IMAGES/traces/",sim,model,"Mars-",sigma,"secs",nyear,"yrs.png")
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
    title=string("IMAGES/traces/",sim,model,"Jupiter-",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
  elseif String(model)=="p4"
    figure(figsize=(8,6))
    for i=1:5
      ax3=subplot(3,2,i)
      for j=1:nwalkers 
      ax3.plot(par_mcmc[j,iburn:nsteps,i+10])
      end
      ax3.set_ylabel(parname[i+10])
    end
    tight_layout()
    title=string("IMAGES/traces/",sim,model,"Mars-",sigma,"secs",nyear,"yrs.png")
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
    title=string("IMAGES/traces/",sim,model,"Jupiter-",sigma,"secs",nyear,"yrs.png")
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
    title=string("IMAGES/traces/",sim,model,"Jupiter-",sigma,"secs",nyear,"yrs.png")
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

function plot_effects(sigma::Real,nyear::Real,sim::String,model::String)
  if String(sim)=="EMB" && isfile(string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  elseif isfile(string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
  file1=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s15yrs.jld2")
  file2=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s16yrs.jld2")
  file3=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s17yrs.jld2")
  file4=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s18yrs.jld2")
  file5=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s19yrs.jld2")
  file6=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s20yrs.jld2")
  file7=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s21yrs.jld2")
  file8=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s22yrs.jld2")
  file9=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s23yrs.jld2")
  file10=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s24yrs.jld2")
  file11=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s25yrs.jld2")
  file12=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s26yrs.jld2")
  file13=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s27yrs.jld2")
  file14=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s28yrs.jld2")
  file15=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s29yrs.jld2")
  file16=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s30yrs.jld2")
  else 
    return  println("MCMC file for ",sim," with ",model," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  m1=jldopen(String(file1),"r")
  m2=jldopen(String(file2),"r")
  m3=jldopen(String(file3),"r")
  m4=jldopen(String(file4),"r")
  m5=jldopen(String(file5),"r")
  m6=jldopen(String(file6),"r")
  m7=jldopen(String(file7),"r")
  m8=jldopen(String(file8),"r")
  m9=jldopen(String(file9),"r")
  m10=jldopen(String(file10),"r")
  m11=jldopen(String(file11),"r")
  m12=jldopen(String(file12),"r")

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