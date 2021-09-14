using PyPlot,Statistics
include("CGS.jl")
rc("font",family="sans-serif")
calc_deg(value) = value * 180/pi
calc_evec1(e,omega) = e* cos(omega-77)
calc_evec2(e,omega) = e* sin(omega-77)
# Basic corner plot for posterior distributions of x vs y parameters
function corner_plot(xvalue,yvalue,nbins,optx,opty,truex,truey)
	meanx=mean(xvalue);sigmax=std(xvalue)
	meany=mean(yvalue);sigmay=std(yvalue)

	fig=figure(figsize=(6,6))
	subplots_adjust(hspace=0.02,wspace=0.02)
	subplot(221)
	ax2=gca()
	ax2.hist(xvalue,bins=nbins,histtype="step",density="true",color="black")
	axvline(meanx-sigmax,color="grey",alpha=0.5,label=L"$1\sigma$ Limit")
	axvline(meanx+sigmax,color="grey",alpha=0.5) 
	axvline(truex,linestyle="--",color="black",label="True Value")
	ax2.minorticks_on()
	ax2.tick_params(which="major",direction="in",length=6,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax2.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax2.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.0)

	subplot(224)
	ax3=gca()
	ax3.hist(yvalue,bins=nbins,histtype="step",density="true",color="black",orientation="horizontal")
	axhline(meany-sigmay,color="grey",alpha=0.5)
	axhline(meany+sigmay,color="grey",alpha=0.5)
	axhline(truey,linestyle="--",color="black")
	ax3.minorticks_on()
	ax3.tick_params(which="major",direction="in",length=6,
	    left="true",right="true",top="false",bottom="false",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax3.tick_params(which="minor",direction="in",length=2,
	    left="true",right="true",top="false",bottom="false",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")

	subplot(223,sharex=ax2,sharey=ax3)
	ax1=gca()
	ax1.hist2d(xvalue,yvalue,bins=nbins,cmin=1)
  axvline(truex,linestyle="--",color="black")
  axhline(truey,linestyle="--",color="black")
	ax1.axis([minimum(xvalue),maximum(xvalue),minimum(yvalue),maximum(yvalue)])
	ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
	ax1.tick_params(which="minor",direction="in",top="true",right="true",length=2)
	show()
end

# Create a corner plot for posterior distributions of moon parameters
function corner_moon(jldmc,nbins) 
  par_mcmc,lprob_mcmc = jldmc["par_mcmc"],jldmc["lprob_mcmc"]
  param = jldmc["param"]
  iburn = jldmc["iburn"]
  nwalkers,nsteps = jldmc["nwalkers"],jldmc["nsteps"]
  x1=vec(par_mcmc[:,iburn:nsteps,16])
  x2=vec(par_mcmc[:,iburn:nsteps,17])
  x3=vec(par_mcmc[:,iburn:nsteps,18])#.*57.2957795
  truex1=0.01
  truex2=0.01
  truex3=2.3122#.*57.2957795

  fig=figure(figsize=(8,8))
  subplots_adjust(hspace=0.05,wspace=0.05)
  subplot(3,3,1)
  ax1=gca()
  ax1.axvline(truex3,linestyle="-",color="black")
  ax1.hist(x3,bins=nbins,histtype="step",density="true",color="black")
  ax1.minorticks_on()
  ax1.tick_params(which="major",direction="in",length=5,
      left="false",right="false",top="true",bottom="true",
      labelbottom="false",labeltop="false",labelleft="false",labelright="false")
  ax1.tick_params(which="minor",direction="in",length=2,
      left="false",right="false",top="true",bottom="true",
      labelbottom="false",labeltop="false",labelleft="false",labelright="false")

  subplot(3,3,4,sharex=ax1)
  ax4=gca()
  ax4.hist2d(x3,x2,bins=nbins,cmin=1)
  ylabel(L"$t_{max} \cos{\phi_0}$")
  ax4.minorticks_on()
  ax4.tick_params(which="major",direction="in",top="true",right="true",length=5,
      labelleft="true",labelbottom="false")
  ax4.tick_params(which="minor",direction="in",top="true",right="true",length=2,
      labelleft="true",labelbottom="false")

  subplot(3,3,5)
  ax5=gca()
  ax5.hist(x2,bins=nbins,histtype="step",density="true",color="black")
  #     ax5.axvline(truex2,linestyle="-",color="black")
  ax5.minorticks_on()
  ax5.tick_params(which="major",direction="in",length=5,
      left="false",right="false",top="true",bottom="true",
      labelbottom="false",labeltop="false",labelleft="false",labelright="false")
  ax5.tick_params(which="minor",direction="in",length=2,
      left="false",right="false",top="true",bottom="true",
      labelbottom="false",labeltop="false",labelleft="false",labelright="false")

  subplot(3,3,7,sharex=ax1)
  ax7=gca()
  ax7.hist2d(x3,x1,bins=nbins,cmin=1)
  xlabel(L"$\Delta \phi [rad]$")
  ylabel(L"$t_{max} \sin{\phi_0}$")
  ax7.minorticks_on()
  ax7.tick_params(which="major",direction="in",top="true",right="true",length=5,
      labelleft="true",labelbottom="true")
  ax7.tick_params(which="minor",direction="in",top="true",right="true",length=2,
      labelleft="true",labelbottom="true")

  subplot(3,3,8,sharex=ax5)
  ax8=gca()
  ax8.hist2d(x2,x1,bins=nbins,cmin=1)
  xlabel(L"$t_{max} \cos{\phi_0}$")
  ax8.minorticks_on()
  ax8.tick_params(which="major",direction="in",top="true",right="true",length=5,
      labelleft="false",labelbottom="true")
  ax8.tick_params(which="minor",direction="in",top="true",right="true",length=2,
      labelleft="false",labelbottom="true")

  subplot(3,3,9)
  ax9=gca()
  ax9.hist(x1,bins=nbins,histtype="step",density="true",color="black")
  xlabel(L"$t_{max} \sin{\phi_0}$")
  #     ax9.axvline(truex1,linestyle="-",color="black")
  #     ylabel(L"$e \cos \varpi $")
  ax9.minorticks_on()
  ax9.tick_params(which="major",direction="in",top="true",right="false",length=5,
      labelleft="false",labelbottom="true")
  ax9.tick_params(which="minor",direction="in",top="true",right="false",length=2,
      labelleft="false",labelbottom="true")
  show()
end
# Create a corner plot for posterior distributions of planet parameters
function corner_planet(sigma,nyear,nbins,pl_name,include_moon::Bool=false) 
    mcfile = string("MCMC/fromEMB/p3_mcmc",sigma,"s",nyear,"yrs.jld2")
    if include_moon
      mcfile = string("MCMC/moon_mcmc",sigma,"s",nyear,"yrs.jld2")
    end
    m = jldopen(String(mcfile),"r")
  par_mcmc= m["par_mcmc"]
  lprob_mcmc = m["lprob_mcmc"]
  nwalkers = m["nwalkers"]
  nsteps = m["nsteps"]
  accept = m["accept"]
  iburn = m["iburn"]
  indepsamples = m["indepsamples"]
  # True values based on "PlanetaryBodyData.pdf" (source?)
  if string(pl_name) == "venus"
    offset = 224.70
    x1=vec(par_mcmc[:,iburn:nsteps,1]).* CGS.MSUN/CGS.MEARTH
    x2=vec(par_mcmc[:,iburn:nsteps,4])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
    x3=vec(par_mcmc[:,iburn:nsteps,5])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
    x4=vec(par_mcmc[:,iburn:nsteps,2]).-offset
    truex1=0.815
    truex2=calc_evec1(0.006,131)
    truex3=calc_evec2(0.006,131)
    truex4=224.7007992.-offset
    lim=0.00076,0.00081
    label=L"Per $- 224.7$ [days]"
  elseif string(pl_name) == "earth"
    offset = 365.25
    x1=vec(par_mcmc[:,iburn:nsteps,6]).* CGS.MSUN/CGS.MEARTH
    x2=vec(par_mcmc[:,iburn:nsteps,9])#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
    x3=vec(par_mcmc[:,iburn:nsteps,10])#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
    x4=vec(par_mcmc[:,iburn:nsteps,7]).-offset
    truex1=1
    truex2=calc_evec1(0.0167,102.4)
    truex3=calc_evec2(0.0167,102.4)
    truex4=365.2564-offset
    lim=0.0064,0.00652
    label=L"Per $- 365.25$ [days]"
  elseif string(pl_name) == "jupiter"
    offset = 0.0
    x1=vec(par_mcmc[:,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH
    x2=vec(par_mcmc[:,iburn:nsteps,14])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    x3=vec(par_mcmc[:,iburn:nsteps,15])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    x4=vec(par_mcmc[:,iburn:nsteps,12])
    truex1=318
    truex2=calc_evec1(0.048,14.75)
    truex3=calc_evec2(0.048,14.75)
    truex4=4332.82012875
    lim=minimum(x4),maximum(x4)
    label="Per [days]"
  end
  println("P_mean= ",mean(x4))
  println(truex2," <- ecosω & esinω -> ",truex3)

	fig=figure(figsize=(9,9))
	subplots_adjust(hspace=0.09,wspace=0.09)
	subplot(4,4,1)
  ax1=gca()
  ax1.axvline(truex4,linestyle="-",color="black")
	ax1.hist(x4,bins=nbins,histtype="step",density="true",color="black")
	ax1.minorticks_on()
	ax1.tick_params(which="major",direction="in",length=5,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax1.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")

  subplot(4,4,5,sharex=ax1)
  ax5=gca()
  ax5.hist2d(x4,x3,bins=nbins,cmin=1)
  ylabel(L"$e \sin \varpi $")
  ax5.minorticks_on()
  ax5.tick_params(which="major",direction="in",top="true",right="true",length=5,
      labelleft="true",labelbottom="false")
  ax5.tick_params(which="minor",direction="in",top="true",right="true",length=2,
      labelleft="true",labelbottom="false")

  subplot(4,4,6)
  ax6=gca()
  ax6.hist(x3,bins=nbins,histtype="step",density="true",color="black")
  ax6.axvline(truex3,linestyle="-",color="black")
  ax6.minorticks_on()
  ax6.tick_params(which="major",direction="in",length=5,
      left="false",right="false",top="true",bottom="true",
      labelbottom="false",labeltop="false",labelleft="false",labelright="false")
  ax6.tick_params(which="minor",direction="in",length=2,
      left="false",right="false",top="true",bottom="true",
      labelbottom="false",labeltop="false",labelleft="false",labelright="false")

  subplot(4,4,9,sharex=ax1)
  ax9=gca()
  ax9.hist2d(x4,x2,bins=nbins,cmin=1)
  ylabel(L"$e \cos \varpi $")
  ax9.minorticks_on()
  ax9.tick_params(which="major",direction="in",top="true",right="true",length=5,
      labelleft="true",labelbottom="false")
  ax9.tick_params(which="minor",direction="in",top="true",right="true",length=2,
      labelleft="true",labelbottom="false")

  subplot(4,4,10,sharex=ax6)
  ax10=gca()
  ax10.hist2d(x3,x2,bins=nbins,cmin=1)
  ax10.minorticks_on()
  ax10.tick_params(which="major",direction="in",top="true",right="true",length=5,
      labelleft="false",labelbottom="false")
  ax10.tick_params(which="minor",direction="in",top="true",right="true",length=2,
      labelleft="false",labelbottom="false")

  subplot(4,4,11)
  ax11=gca()
  ax11.hist(x2,bins=nbins,histtype="step",density="true",color="black")
  ax11.axvline(truex2,linestyle="-",color="black")
  ax11.minorticks_on()
  ax11.tick_params(which="major",direction="in",length=5,
      left="false",right="false",top="true",bottom="true",
      labelbottom="false",labeltop="false",labelleft="false",labelright="false")
  ax11.tick_params(which="minor",direction="in",length=2,
      left="false",right="false",top="true",bottom="true",
      labelbottom="false",labeltop="false",labelleft="false",labelright="false")    

  subplot(4,4,13,sharex=ax1)
  ax13=gca()
  ax13.hist2d(x4,x1,bins=nbins,cmin=1)
  # xlim(lim)
  xlabel(label)
  ylabel(L"Mass [$M_{Earth}$]")
  ax13.minorticks_on()
  ax13.tick_params(which="major",direction="in",top="true",right="true",length=5,
      labelbottom="true",labelleft="true")
  ax13.tick_params(which="minor",direction="in",top="true",right="true",length=2,
      labelbottom="false",labelleft="true")

  subplot(4,4,14,sharex=ax6)
  ax14=gca()
  ax14.hist2d(x3,x1,bins=nbins,cmin=1)
  xlabel(L"$e \sin \varpi $")
  ax14.minorticks_on()
  ax14.tick_params(which="major",direction="in",top="true",right="true",length=5,
      labelleft="false",labelbottom="true")
  ax14.tick_params(which="minor",direction="in",top="true",right="true",length=2,
      labelleft="false",labelbottom="false")

  subplot(4,4,15,sharex=ax11)
  ax15=gca()
  ax15.hist2d(x2,x1,bins=nbins,cmin=1)
  xlabel(L"$e \cos \varpi $")
  ax15.minorticks_on()
  ax15.tick_params(which="major",direction="in",top="true",right="true",length=5,
      labelleft="false",labelbottom="true")
  ax15.tick_params(which="minor",direction="in",top="true",right="true",length=2,
      labelleft="false",labelbottom="false")

  subplot(4,4,16)
  ax16=gca()
  ax16.hist(x1,bins=nbins,histtype="step",density="true",color="black")
  ax16.axvline(truex1,linestyle="-",color="black")
  xlabel(L"Mass [$M_{Earth}$]")
  ax16.minorticks_on()
  ax16.tick_params(which="major",direction="in",length=5,
      left="false",right="false",top="true",bottom="true",
      labelbottom="true",labeltop="false",labelleft="false",labelright="false")
  ax16.tick_params(which="minor",direction="in",length=2,
      left="false",right="false",top="true",bottom="true",
      labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	show()
end