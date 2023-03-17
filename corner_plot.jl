using JLD2,PyPlot,Statistics,Distributions,LinearAlgebra
using3D() # Needed to create a 3D subplot
#include("CGS.jl")
include("MCMC.jl")
rc("font",family="sans-serif")
rc("lines",linewidth=1.5)
calc_deg(value)=value * 180/pi
calc_evec1(e,omega)=e* cos(omega-77)
calc_evec2(e,omega)=e* sin(omega-77)
calc_tmax(a_p,a_s,m_p,m_s,P_p)=(a_s*m_s*P_p) / (2*pi*a_p*(m_s+m_p))
# Basic corner plot for posterior distributions of 2 parameters
function corner(x1,x2,nbins)
  fig=figure(figsize=(5,5),dpi=150)
  subplots_adjust(hspace=0.05,wspace=0.05)
  ax2=subplot(221)
  ax2.hist(x1,bins=nbins,histtype="step",color="black")
  ax2.minorticks_on()
  ax2.tick_params(which="both",direction="out",
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax3=subplot(224)
  ax3.hist(x2,bins=nbins,histtype="step",color="black",orientation="horizontal")
  ax3.minorticks_on()
  ax3.tick_params(which="both",direction="out",
      left=true,right=true,top=false,bottom=false,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax1=subplot(223,sharex=ax2,sharey=ax3)
  ax1.hist2d(x1,x2,bins=nbins,cmin=1)
  ax1.axis([minimum(x1),maximum(x1),minimum(x2),maximum(x2)])
  ax1.tick_params(which="both",direction="out",top=true,right=true)
  return fig
end
# Corner plot for posterior distributions of 2 parameters, compared to true values
function corner(x1,x2,truex1,truex2,nbins)
	fig=figure(figsize=(4,4),dpi=150)
	# subplots_adjust(hspace=0.05,wspace=0.05)
  # ax2=fig.add_axes(221)
  # text(0.5,0.5,"test")
  ax1=fig.add_subplot(223)
  ax1.hist2d(x1,x2,bins=nbins,cmin=1,cmax=100000)
  ax1.axis([minimum(x1),maximum(x1),minimum(x2),maximum(x2)])
  xlabel(L"$t_{max}$ [min], TTV Amplitude on Planet 2")
  ylabel(L"$\Delta\phi$ [rad], $Per_2$ Phase Offset")
  ax1.tick_params(which="both",direction="out",top=true,right=true)

	ax2=fig.add_subplot(221,sharex=ax1)
	h1=ax2.hist(x1,bins=nbins,histtype="step",color="black")
  ax2.axvspan(quantile(x1,0.1587),quantile(x1,0.8413),color="limegreen",alpha=0.35)
	axvline(truex1,linestyle="--",color="black",label="True Value")
	ax2.minorticks_on()
	ax2.tick_params(which="both",direction="out",
	    left=false,right=false,top=false,bottom=true,
	    labelbottom=true,labeltop=false,labelleft=false,labelright=false)

	ax3=fig.add_subplot(224,sharey=ax1)
	h2=ax3.hist(x2,bins=nbins,histtype="step",color="black",orientation="horizontal")
  ax3.axhspan(quantile(x2,0.1587),quantile(x2,0.8413),color="limegreen",alpha=0.35)
	axhline(truex2,linestyle="--",color="black")
	ax3.minorticks_on()
	ax3.tick_params(which="both",direction="out",
	    left=false,right=true,top=false,bottom=false,
	    labelbottom=false,labeltop=false,labelleft=false,labelright=true)
  return fig
end
# Corner plot for posterior distributions of 3 parameters
function corner(x1,x2,x3,nbins)
  fig=figure(figsize=(5,5),dpi=150)
  subplots_adjust(hspace=0.05,wspace=0.05)
  ax1=subplot(3,3,1)
  ax1.hist(x3,bins=nbins,histtype="step",color="black")
  ax1.minorticks_on()
  ax1.tick_params(which="major",direction="out",length=5,
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)
  ax1.tick_params(which="minor",direction="out",length=2,
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax4=subplot(3,3,4,sharex=ax1)
  ax4.hist2d(x3,x2,bins=nbins,cmin=1)
  ylabel(L"$t_{max} \cos{\phi_0}$")
  ax4.minorticks_on()
  ax4.tick_params(which="major",direction="out",top=true,right=true,length=5,
      labelleft=true,labelbottom=false)
  ax4.tick_params(which="minor",direction="out",top=true,right=true,length=2,
      labelleft=true,labelbottom=false)

  ax5=subplot(3,3,5)
  ax5.hist(x2,bins=nbins,histtype="step",color="black")
  ax5.minorticks_on()
  ax5.tick_params(which="major",direction="out",length=5,
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)
  ax5.tick_params(which="minor",direction="out",length=2,
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax7=subplot(3,3,7,sharex=ax1)
  ax7.hist2d(x3,x1,bins=nbins,cmin=1)
  xlabel(L"$\Delta \phi [rad]$")
  ylabel(L"$t_{max} \sin{\phi_0}$")
  ax7.minorticks_on()
  ax7.tick_params(which="major",direction="out",top=true,right=true,length=5,
      labelleft=true,labelbottom=true)
  ax7.tick_params(which="minor",direction="out",top=true,right=true,length=2,
      labelleft=true,labelbottom=true)

  ax8=subplot(3,3,8,sharex=ax5)
  ax8.hist2d(x2,x1,bins=nbins,cmin=1)
  xlabel(L"$t_{max} \cos{\phi_0}$")
  ax8.minorticks_on()
  ax8.tick_params(which="major",direction="out",top=true,right=true,length=5,
      labelleft=false,labelbottom=true)
  ax8.tick_params(which="minor",direction="out",top=true,right=true,length=2,
      labelleft=false,labelbottom=true)

  ax9=subplot(3,3,9)
  ax9.hist(x1,bins=nbins,histtype="step",color="black")
  xlabel(L"$t_{max} \sin{\phi_0}$")
  ax9.minorticks_on()
  ax9.tick_params(which="major",direction="out",top=true,left=false,right=false,length=5,
      labelleft=false,labelbottom=true)
  ax9.tick_params(which="minor",direction="out",top=true,left=false,right=false,length=2,
      labelleft=false,labelbottom=true)
  tight_layout()
end
# Corner plot for posterior distributions of 3 parameters, compared to true values
function corner(x1,x2,x3,truex1,truex2,truex3,nbins)
  fig=figure(figsize=(6,6),dpi=150)
  subplots_adjust(hspace=0.05,wspace=0.05)
  ax1=subplot(3,3,1)
  ax1.hist(x3,bins=nbins,histtype="step",color="black")
  ax1.axvline(truex3,linestyle="-",color="black")
  ax1.axvspan(quantile(x3,0.1587),quantile(x3,0.8413),color="limegreen",alpha=0.5)
  ax1.minorticks_on()
  ax1.tick_params(which="both",direction="out",
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax4=subplot(3,3,4,sharex=ax1)
  ax4.hist2d(x3,x2,bins=nbins,cmin=1)
  ylabel(L"$t_{max} \cos{\phi_0}$")
  ax4.minorticks_on()
  ax4.tick_params(which="both",direction="out",top=true,right=true,
      labelleft=true,labelbottom=false)

  ax5=subplot(3,3,5)
  ax5.hist(x2,bins=nbins,histtype="step",color="black")
  ax5.axvline(truex2,linestyle="-",color="black")
  ax5.axvspan(quantile(x2,0.1587),quantile(x2,0.8413),color="limegreen",alpha=0.5)
  ax5.minorticks_on()
  ax5.tick_params(which="both",direction="out",
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax7=subplot(3,3,7,sharex=ax1)
  ax7.hist2d(x3,x1,bins=nbins,cmin=1)
  xlabel(L"$\Delta \phi [rad]$")
  ylabel(L"$t_{max} \sin{\phi_0}$")
  ax7.minorticks_on()
  ax7.tick_params(which="both",direction="out",top=true,right=true,
      labelleft=true,labelbottom=true)

  ax8=subplot(3,3,8,sharex=ax5)
  ax8.hist2d(x2,x1,bins=nbins,cmin=1)
  xlabel(L"$t_{max} \cos{\phi_0}$")
  ax8.minorticks_on()
  ax8.tick_params(which="both",direction="out",top=true,right=true,
      labelleft=false,labelbottom=true)

  ax9=subplot(3,3,9)
  ax9.hist(x1,bins=nbins,histtype="step",color="black")
  ax9.axvline(truex1,linestyle="-",color="black")
  ax9.axvspan(quantile(x1,0.1587),quantile(x1,0.8413),color="limegreen",alpha=0.5)
  xlabel(L"$t_{max} \sin{\phi_0}$")
  #     ylabel(L"$e \cos \varpi $")
  ax9.minorticks_on()
  ax9.tick_params(which="both",direction="out",top=true,left=false,right=false,
      labelleft=false,labelbottom=true)
  tight_layout()
  return fig
end
# Corner plot for posterior distributions of 4 parameters, compared to true values
function corner(x1,x2,x3,x4,truex1,truex2,truex3,truex4,nbins,lim,label)
  fig=figure(figsize=(6.5,6.5),dpi=150)
  subplots_adjust(hspace=0.05,wspace=0.05)
  ax1=fig.add_subplot(4,4,1)
  ax1.axvline(truex4,linestyle="--",color="black")
  ax1.hist(x4,bins=nbins,histtype="step",color="black")
  ax1.axvspan(quantile(x4,0.1587),quantile(x4,0.8413),color="limegreen",alpha=0.35)
  ax1.minorticks_on()
  ax1.tick_params(which="both",direction="in",
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax5=fig.add_subplot(4,4,5,sharex=ax1)
  h21=ax5.hist2d(x4,x3,bins=nbins,cmin=1,cmax=100000)
  ylabel(L"$e \sin \varpi $")
  ax5.minorticks_on()
  ax5.tick_params(which="both",direction="out",top=true,right=false,
      labelleft=true,labelbottom=false)

  ax6=fig.add_subplot(4,4,6)
  ax6.hist(x3,bins=nbins,histtype="step",color="black")
  ax6.axvline(truex3,linestyle="--",color="black")
  ax6.axvspan(quantile(x3,0.1587),quantile(x3,0.8413),color="limegreen",alpha=0.35)
  ax6.minorticks_on()
  ax6.tick_params(which="both",direction="in",
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax9=fig.add_subplot(4,4,9,sharex=ax1)
  h22=ax9.hist2d(x4,x2,bins=nbins,cmin=1,cmax=100000)
  ylabel(L"$e \cos \varpi $")
  ax9.minorticks_on()
  ax9.tick_params(which="both",direction="out",top=true,right=false,
      labelleft=true,labelbottom=false)

  ax10=fig.add_subplot(4,4,10,sharex=ax6)
  h23=ax10.hist2d(x3,x2,bins=nbins,cmin=1,cmax=100000)
  ax10.minorticks_on()
  ax10.tick_params(which="both",direction="out",top=true,right=false,
      labelleft=false,labelbottom=false)

  ax11=fig.add_subplot(4,4,11)
  ax11.hist(x2,bins=nbins,histtype="step",color="black")
  ax11.axvline(truex2,linestyle="--",color="black")
  ax11.axvspan(quantile(x2,0.1587),quantile(x2,0.8413),color="limegreen",alpha=0.35)
  ax11.minorticks_on()
  ax11.tick_params(which="both",direction="in",
      left=false,right=false,top=true,bottom=true,
      labelbottom=false,labeltop=false,labelleft=false,labelright=false) 

  ax13=fig.add_subplot(4,4,13,sharex=ax1)
  h24=ax13.hist2d(x4,x1,bins=nbins,cmin=1,cmax=100000)
  xlim(lim)
  xlabel(label)
  ylabel(L"Mass [$M_{\oplus}$]")
  ax13.minorticks_on()
  ax13.tick_params(which="both",direction="out",top=true,right=false,
      labelbottom=true,labelleft=true)

  ax14=fig.add_subplot(4,4,14,sharex=ax6)
  h25=ax14.hist2d(x3,x1,bins=nbins,cmin=1,cmax=100000)
  xlabel(L"$e \sin \varpi $")
  ax14.minorticks_on()
  ax14.tick_params(which="both",direction="out",top=true,right=false,
      labelleft=false,labelbottom=true)

  ax15=fig.add_subplot(4,4,15,sharex=ax11)
  h26=ax15.hist2d(x2,x1,bins=nbins,cmin=1,cmax=100000)
  xlabel(L"$e \cos \varpi $")
  ax15.minorticks_on()
  ax15.tick_params(which="both",direction="out",top=true,right=false,
      labelleft=false,labelbottom=true)

  ax16=fig.add_subplot(4,4,16)
  ax16.hist(x1,bins=nbins,histtype="step",color="black")
  ax16.axvline(truex1,linestyle="--",color="black")
  ax16.axvspan(quantile(x1,0.1587),quantile(x1,0.8413),color="limegreen",alpha=0.35)
  xlabel(L"Mass [$M_{\oplus}$]")
  ax16.minorticks_on()
  ax16.tick_params(which="both",direction="in",
      left=false,right=false,top=true,bottom=true,
      labelbottom=true,labeltop=false,labelleft=false,labelright=false)
  tight_layout()
 return fig
end
# Create a corner plot for significant posterior distributions of planet parameters
function corner_hist(sigma,nyear,grid_type_nplanet,case_num,nbins,include_moon::Bool=false) 
  if case_num==1  && isfile(string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  elseif case_num==2 && isfile(string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  else
    return  println("MCMC file for case: ",case_num," with ",grid_type_nplanet," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  m=jldopen(String(mcfile),"r")
  par_mcmc= m["par_mcmc"]
  lprob_mcmc=m["lprob_mcmc"]
  nwalkers=m["nwalkers"]
  nsteps=m["nsteps"]
  iburn=m["iburn"]
  indepsamples=m["indepsamples"]
  # True values based on "PlanetaryBodyData.pdf" (source?)
  offset=224.70
  m1=vec(par_mcmc[:,iburn:nsteps,1]).* CGS.MSUN/CGS.MEARTH
  ec1=vec(par_mcmc[:,iburn:nsteps,4])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
  es1=vec(par_mcmc[:,iburn:nsteps,5])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
  p1=vec(par_mcmc[:,iburn:nsteps,2]).-offset
  e1=sqrt.(vec(par_mcmc[:,iburn:nsteps,4]).^2 .+ vec(par_mcmc[:,iburn:nsteps,5]).^2)
  truem1=0.815
  trueec1=calc_evec1(0.00677323,131.53298)
  truees1=calc_evec2(0.00677323,131.53298)
  truep1=224.7007992.-offset
  truee1=0.00677323
  lim=0.00076,0.00081
  offset=365.25
  m2=vec(par_mcmc[:,iburn:nsteps,6]).* CGS.MSUN/CGS.MEARTH
  ec2=vec(par_mcmc[:,iburn:nsteps,9])
  es2=vec(par_mcmc[:,iburn:nsteps,10])#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
  p2=vec(par_mcmc[:,iburn:nsteps,7]).-offset
  e2=sqrt.(vec(par_mcmc[:,iburn:nsteps,9]).^2 .+ vec(par_mcmc[:,iburn:nsteps,10]).^2)
  truem2=1
  trueec2=calc_evec1(0.01671022,102.94719)
  truees2=calc_evec2(0.01671022,102.94719)
  truep2=365.2564-offset #365.256355
  truee2=0.01671022
  lim=0.0064,0.00652
  corner(m1,m2,truem1,truem2,nbins)
  ylabel(L"Mass of Earth [$M_{Earth}$]")
  xlabel(L"Mass of Venus [$M_{Earth}$]")
  tight_layout()
  title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"masses",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  corner(e1,e2,truee1,truee2,nbins)
  ylabel("Eccentricity of Earth")
  xlabel("Eccentricity of Venus")
  tight_layout()
  title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"eccs",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  corner(ec1,ec2,trueec1,trueec2,nbins)
  ylabel(L"$e \cos \varpi$ for Earth")
  xlabel(L"$e \cos \varpi$ for Venus")
  tight_layout()
  title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"ecos",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  corner(es1,es2,truees1,truees2,nbins)
  ylabel(L"$e \sin \varpi$ for Earth")
  xlabel(L"$e \sin \varpi$ for Venus")
  tight_layout()
  title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"esin",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  if grid_type_nplanet=="p4"
    m3=vec(par_mcmc[:,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH
    ec3=vec(par_mcmc[:,iburn:nsteps,14])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    es3=vec(par_mcmc[:,iburn:nsteps,15])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    p3=vec(par_mcmc[:,iburn:nsteps,12])
    e3=sqrt.(vec(par_mcmc[:,iburn:nsteps,14]).^2 .+ vec(par_mcmc[:,iburn:nsteps,15]).^2)
    truem3=0.1074
    trueec3=calc_evec1(0.09341233,336.04084)
    truees3=calc_evec2(0.09341233,336.04084)
    truep3=686.9795859
    truee3=0.09341233
    corner(m3,e3,truem3,truee3,nbins)
    xlabel(L"Mass of Mars [$M_{Earth}$]")
    ylabel("Eccentricity of Mars")
    tight_layout()
    title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"Vmecc",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
    m4=vec(par_mcmc[:,iburn:nsteps,16]).* CGS.MSUN/CGS.MEARTH
    ec4=vec(par_mcmc[:,iburn:nsteps,19])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    es4=vec(par_mcmc[:,iburn:nsteps,20])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    p4=vec(par_mcmc[:,iburn:nsteps,17])
    e4=sqrt.(vec(par_mcmc[:,iburn:nsteps,19]).^2 .+ vec(par_mcmc[:,iburn:nsteps,20]).^2)
    truem4=317.8
    trueec4=calc_evec1(0.04839266,14.75385)
    truees4=calc_evec2(0.04839266,14.75385)
    truep4=4332.82012875
    truee4=0.04839266
    corner(m4,e4,truem4,truee4,nbins)
    xlabel(L"Mass of Jupiter [$M_{Earth}$]")
    ylabel("Eccentricity of Jupiter")
    tight_layout()
    title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"Jmecc",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
  else 
    m3=vec(par_mcmc[:,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH
    ec3=vec(par_mcmc[:,iburn:nsteps,14])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    es3=vec(par_mcmc[:,iburn:nsteps,15])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    p3=vec(par_mcmc[:,iburn:nsteps,12])
    e3=sqrt.(vec(par_mcmc[:,iburn:nsteps,14]).^2 .+ vec(par_mcmc[:,iburn:nsteps,15]).^2)
    truem3=317.8
    trueec3=calc_evec1(0.04839266,14.75385)
    truees3=calc_evec2(0.04839266,14.75385)
    truep3=4332.82012875
    truee3=0.04839266
    corner(m3,e3,truem3,truee3,nbins)
    xlabel(L"Mass of Jupiter [$M_{Earth}$]")
    ylabel("Eccentricity of Jupiter")
    tight_layout()
    title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"mecc",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
  end
  if include_moon
    tmax=sqrt.(vec(par_mcmc[:,iburn:nsteps,16]).^2 .+ vec(par_mcmc[:,iburn:nsteps,17]).^2)
    x1=vec(par_mcmc[:,iburn:nsteps,16])
    x2=vec(par_mcmc[:,iburn:nsteps,17])
    x3=vec(par_mcmc[:,iburn:nsteps,18])#.*57.2957795
    truetmax=calc_tmax(CGS.AU,CGS.AMOON*CGS.AU,CGS.MEARTH,CGS.MMOON,365.256355) #0.0018
    truex2=0.01
    truex3=2.31586#.*57.2957795
    title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"Moon2-",sigma,"secs",nyear,"yrs.png")
    corner(tmax,x3,truetmax,truex3,nbins)
    xlabel(L"$t_{max}$ [days]")
    ylabel(L"$\Delta \phi$ [rad]")
    tight_layout()
    savefig(title)
    clf()
  end
  # show()
end
# Create a corner plot for posterior distributions of planet parameters
function corner_plot(sigma,nyear,grid_type_nplanet,case_num,nbins,include_moon::Bool=false) 
  if case_num==1  #&& isfile(string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  elseif case_num==2 #&& isfile(string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  else
    return  println("MCMC file for case",case_num," with ",grid_type_nplanet," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  m=jldopen(String(mcfile),"r")
  par_mcmc= m["par_mcmc"]
  lprob_mcmc=m["lprob_mcmc"]
  nwalkers=m["nwalkers"]
  nsteps=m["nsteps"]
  iburn=m["iburn"]
  indepsamples=m["indepsamples"]
  model_BIC,sigsys,sigtot,chi=mc_vals(sigma,nyear,grid_type_nplanet,case_num,include_moon)
  # println("Loaded",mcfile,".")
  if grid_type_nplanet=="p3" || grid_type_nplanet=="widep3"
    model=L"$\mathcal{H}_{PPPP}$"
  elseif grid_type_nplanet=="p4" || grid_type_nplanet=="widep4"
    model=L"$\mathcal{H}_{PPPP}$"
  elseif grid_type_nplanet=="p3moon" || grid_type_nplanet=="widep3moon"
    model=L"$\mathcal{H}_{PPsP}$"
  elseif grid_type_nplanet=="p3moonp4" || grid_type_nplanet=="widep3moonp4"
    model=L"$\mathcal{H}_{PPsPP}$"
  end
  # True values based on "PlanetaryBodyData.pdf" (source?)
  offset=224.70
  m1=vec(par_mcmc[:,iburn:nsteps,1]).* CGS.MSUN/CGS.MEARTH
  ec1=vec(par_mcmc[:,iburn:nsteps,4])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
  es1=vec(par_mcmc[:,iburn:nsteps,5])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
  p1=vec(par_mcmc[:,iburn:nsteps,2]).-offset
  truem1=0.815
  trueec1=calc_evec1(0.00677323,131.53298)
  truees1=calc_evec2(0.00677323,131.53298)
  truep1=224.7007992.-offset
  lim=0.00076,0.00081
  label=L"Per$_1 - 224.7$ [days]"
  title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Venus-",sigma,"secs",nyear,"yrs.png")
  fig1=corner(m1,ec1,es1,p1,truem1,trueec1,truees1,truep1,nbins,lim,label)
  fig1.suptitle(string(model," Posteriors for Planet 1"))
  fig1.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sig_tot," sec",'\n',"BIC= ",model_BIC,'\n',L"$\chi^2 =$",chi))
  savefig(title)
  clf()
  offset=365.25
  m2=vec(par_mcmc[:,iburn:nsteps,6]).* CGS.MSUN/CGS.MEARTH;
  ec2=vec(par_mcmc[:,iburn:nsteps,9]);#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
  es2=vec(par_mcmc[:,iburn:nsteps,10]);#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
  p2=vec(par_mcmc[:,iburn:nsteps,7]).-offset;
  truem2=1
  trueec2=calc_evec1(0.01671022,102.94719)
  truees2=calc_evec2(0.01671022,102.94719)
  truep2=365.2564-offset #365.256355
  lim=minimum(p2),maximum(p2)#0.0064,0.00652
  label=L"Per$_2 - 365.25$ [days]"
  title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Earth-",sigma,"secs",nyear,"yrs.png")
  fig2=corner(m2,ec2,es2,p2,truem2,trueec2,truees2,truep2,nbins,lim,label)
  fig2.suptitle(string(model," Posteriors for Planet 2"))
  fig2.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sig_tot," sec",'\n',"BIC= ",model_BIC,'\n',L"$\chi^2 =$",chi))
  savefig(title)
  clf()
  if grid_type_nplanet=="p4" || grid_type_nplanet=="p3moonp4" || grid_type_nplanet=="widep4"
    m3=vec(par_mcmc[:,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH
    ec3=vec(par_mcmc[:,iburn:nsteps,14])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    es3=vec(par_mcmc[:,iburn:nsteps,15])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    p3=vec(par_mcmc[:,iburn:nsteps,12])
    truem3=0.1074
    trueec3=calc_evec1(0.09341233,336.04084)
    truees3=calc_evec2(0.09341233,336.04084)
    truep3=686.9795859
    lim=minimum(p3),maximum(p3)
    label=L"Per$_4$ [days]"
    title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Mars-",sigma,"secs",nyear,"yrs.png")
    fig3=corner(m3,ec3,es3,p3,truem3,trueec3,truees3,truep3,nbins,lim,label)
    fig3.suptitle(string(model," Posteriors for Planet 4"))
    fig3.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sig_tot," sec",'\n',"BIC= ",model_BIC,'\n',L"$\chi^2 =$",chi))
    savefig(title)
    clf()
    #  println("Mars= ",minimum(m3)," ",minimum(p3)," ",minimum(ec3)," ",minimum(es3))
    #  println("Mars= ",maximum(m3)," ",maximum(p3)," ",maximum(ec3)," ",maximum(es3))
    m4=vec(par_mcmc[:,iburn:nsteps,16]).* CGS.MSUN/CGS.MEARTH
    ec4=vec(par_mcmc[:,iburn:nsteps,19])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    es4=vec(par_mcmc[:,iburn:nsteps,20])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    p4=vec(par_mcmc[:,iburn:nsteps,17])
    truem4=317.8
    trueec4=calc_evec1(0.04839266,14.75385)
    truees4=calc_evec2(0.04839266,14.75385)
    truep4=4332.82012875
    lim=minimum(p4),maximum(p4)
    label=L"Per$_3$ [days]"
    title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Jupiter-",sigma,"secs",nyear,"yrs.png")
     # println("Jupiter= ",minimum(m4)," ",minimum(p4)," ",minimum(ec4)," ",minimum(es4))
     # println("Jupiter= ",maximum(m4)," ",maximum(p4)," ",maximum(ec4)," ",maximum(es4))
    fig4=corner(m4,ec4,es4,p4,truem4,trueec4,truees4,truep4,nbins,lim,label)
    fig4.suptitle(string(model," Posteriors for Planet 3"))
    fig4.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sig_tot," sec",'\n',"BIC= ",model_BIC,'\n',L"$\chi^2 =$",chi))
    savefig(title)
    clf()
  else
    m3=vec(par_mcmc[:,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH
    ec3=vec(par_mcmc[:,iburn:nsteps,14])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    es3=vec(par_mcmc[:,iburn:nsteps,15])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    p3=vec(par_mcmc[:,iburn:nsteps,12])
    truem3=317.8
    trueec3=calc_evec1(0.04839266,14.75385)
    truees3=calc_evec2(0.04839266,14.75385)
    truep3=4332.82012875
    lim=minimum(p3),maximum(p3)
    label=L"Per$_3$ [days]"
    title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Jupiter-",sigma,"secs",nyear,"yrs.png")
     # println("Jupiter= ",minimum(m3)," ",minimum(p3)," ",minimum(ec3)," ",minimum(es3))
     # println("Jupiter= ",maximum(m3)," ",maximum(p3)," ",maximum(ec3)," ",maximum(es3))
    fig3=corner(m3,ec3,es3,p3,truem3,trueec3,truees3,truep3,nbins,lim,label)
    fig3.suptitle(string(model," Posteriors for Planet 3"))
    fig3.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sig_tot," sec",'\n',"BIC= ",model_BIC,'\n',L"$\chi^2 =$",chi))
    savefig(title)
    clf()
  end
  if include_moon && grid_type_nplanet=="p3moon"
     tmax=sqrt.(vec(par_mcmc[:,iburn:nsteps,16]).^2 .+ vec(par_mcmc[:,iburn:nsteps,17]).^2)
     x1=vec(par_mcmc[:,iburn:nsteps,16])
     x2=vec(par_mcmc[:,iburn:nsteps,17])
     x3=vec(par_mcmc[:,iburn:nsteps,18])#.*57.2957795
     truetmax=calc_tmax(CGS.AU,CGS.AMOON*CGS.AU,CGS.MEARTH,CGS.MMOON,365.256355) #0.0018
     truex2=0.01
     truex3=2.31586#.*57.2957795
     title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Moon2-",sigma,"secs",nyear,"yrs.png")
     fig5=corner(tmax,x3,truetmax,truex3,nbins)
     fig5.suptitle(string(model," Posteriors for Satellite"))
     fig5.text(0.575,0.725,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sig_tot," sec",'\n',"BIC= ",model_BIC,'\n',L"$\chi^2 =$",chi))
  #   corner(x1,x2,x3,nbins)
       savefig(title)
       clf()
       elseif include_moon 
     tmax=sqrt.(vec(par_mcmc[:,iburn:nsteps,21]).^2 .+ vec(par_mcmc[:,iburn:nsteps,22]).^2) .* 24*60
     x1=vec(par_mcmc[:,iburn:nsteps,21])
     x2=vec(par_mcmc[:,iburn:nsteps,22])
     x3=vec(par_mcmc[:,iburn:nsteps,23])#.*57.2957795
     truetmax=calc_tmax(CGS.AU,CGS.AMOON*CGS.AU,CGS.MEARTH,CGS.MMOON,365.256355) .* 24*60#0.0018
     truex2=0.01
     truex3=2.31586#.*57.2957795
     title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Moon2-",sigma,"secs",nyear,"yrs.png")
     fig5=corner(tmax,x3,truetmax,truex3,nbins)
     fig5.suptitle(string(model," Posteriors for Satellite"))
     fig5.text(0.575,0.725,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sig_tot," sec",'\n',"BIC= ",model_BIC,'\n',L"$\chi^2 =$",chi))
  #   corner(x1,x2,x3,nbins)
       savefig(title)
       clf()
       end
  # close()
end