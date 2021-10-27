using PyPlot,Statistics,Distributions,Optim,JLD2
rc("font",family="sans-serif")
include("histogram.jl")
avg(x,y) = (x + y)/2
gaussian(x,mu,sig)=exp.(-((x .- mu).^2) ./ (2 * sig^.2))
# Create a plot of likelihood/probability of values in grid
# function plot_prob()
	
# end
function plot_likelihood(sigma,nyear,sim,fitmodel,mcmodel,nbins,include_moon::Bool=false)
	if String(sim)=="EMB" && isfile(string("MCMC/fromEMB/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")) 
    mcfile = string("MCMC/fromEMB/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile = string("FITS/fromEMB/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
  elseif isfile(string("MCMC/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")) 
    mcfile = string("MCMC/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile = string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
    println("Loading...",mcfile," and ",fitfile)
  else 
    return  println("MCMC or FITS file for ",sim," with ",mcmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
 
	wide = jldopen("FITS/wide2_fit30.0s30.0yrs.jld2","r")
	grid_wide = wide["p3"]./365.25
	lprob_wide = exp.((wide["lprob_p3"] .-maximum(wide["lprob_p3"])))
	m = jldopen(String(mcfile),"r")
	f = jldopen(String(fitfile),"r")

	xgrid = f["p3"]./365.25
	xprob = exp.((f["lprob_p3"] .-maximum(f["lprob_p3"])))
	grid = [grid_wide[1:155];xgrid;grid_wide[198:end]]
	prob = [lprob_wide[1:155];xprob;lprob_wide[198:end]]
	# println("Simulated with σ= ",sigma," second noise, for ",nyear," observation years.")
	truex = 11.862615
	ncol = 12
	bfvalue = f["best_p3"][ncol] /365.25	
	param = vec(m["par_mcmc"][:,m["iburn"]:end,ncol])./365.25
	xbin,xhist,xbin_square,hist_square=histogram(param,nbins)
	label="Planet Period Search Grid [years]"
	lim = minimum(xgrid),maximum(xgrid)
	fig = plt.figure()
  ax1 = gca()
	ax1.axvline(truex,linestyle="--",color="black",label=string(truex))
	ax1.text(12,1,"Jupiter")
	# ax1.plot(grid_wide,lprob_wide,linestyle="-",color="grey")
	ax1.plot(grid,prob,color="black",alpha=0.75)
	# xlim(minimum(grid_wide),14)
	# ax1.legend(loc="upper left")
	xlabel(label)
	ylabel("Probability")
	ax1.minorticks_on()
	ax1.tick_params(which="both",direction="in")
  # savefig(title)
  # clf()
  if fitmodel=="p4"
	  xgrid4 = f["p4"]./365.25
		xprob4 = exp.((f["lprob_p4"] .-maximum(f["lprob_p4"])))
		# grid = [grid_wide[1:55];xgrid;grid_wide[60:end]]
		# prob = [lprob_wide[1:55];xprob;lprob_wide[60:end]]
		bfvalue4 = f["best_p4"][ncol] /365.25	
		truex4 = 1.88
		param4 = vec(m["par_mcmc"][:,m["iburn"]:end,ncol])./365.25
		xbin4,xhist4,xbin_square4,hist_square4=histogram(param4,nbins)
		ax1.plot(xgrid,xprob,color="firebrick") 
		ax1.axvline(truex4,linestyle="--",color="black",label=string(truex4))
		ax1.plot(xgrid4,xprob4,color="darkcyan")
		ax1.text(2,1,"Mars")
		ax1.set_xlim(1,14.3)
	  title=string("IMAGES/likelihoods/",sim,mcmodel,"Mars-",sigma,"secs",nyear,"yrs.png")
	else
		ax1.plot(grid,prob,color="firebrick") 
		# Inset zoom of finer period grid 
	  ax2 = fig.add_axes([0.2,0.4,0.4,0.5])
	  # ax2.axvline(bfvalue,linestyle="--",color="black")
	  ax2.plot(xgrid,xprob,color="firebrick")
	  ax2.plot(xbin_square,hist_square./maximum(hist_square),color="firebrick",linewidth=2,alpha=0.5)
	  xlim(lim)
		# ax2.minorticks_on()
	  ax2.tick_params(which="both",direction="in",left="true",labelleft="true")
	  title=string("IMAGES/likelihoods/",sim,mcmodel,"Jupiter-",sigma,"secs",nyear,"yrs.png")
	  tight_layout()
  end
	if include_moon
		xgrid = range(f["dpin"],stop=f["dpout"],length=f["ndp"])
		xprob = exp.((f["lprob_dp"] .-maximum(f["lprob_dp"])))
		grid_wide = range(wide["dpin"],stop=wide["dpout"],length=wide["ndp"])
		lprob_wide = exp.((wide["lprob_dp"] .-maximum(wide["lprob_dp"])))
		grid = [grid_wide[1:60];xgrid;grid_wide[75:end]]
		prob = [lprob_wide[1:60];xprob;lprob_wide[75:end]]
		truex = 2.31586
		ncol = 18
		bfvalue = f["best_dp"][ncol]
		param = vec(m["par_mcmc"][:,m["iburn"]:end,ncol])
		xbin,xhist,xbin_square,hist_square=histogram(param,nbins)
		label="Satellite Phase Offset Search Grid [radians]"
		lim = f["dpin"],f["dpout"]
		fig = plt.figure()
	  ax3 = gca()
		ax3.axvline(truex,linestyle="--",color="black",label=string(truex))
		ax3.text(2.4,1,"Moon")
		ax3.plot(grid,prob,color="purple") 
		# tight_layout()
		# ax3.plot(grid_wide,lprob_wide,color="red") 
		xlabel(label)
		ylabel("Probability")
		ax3.minorticks_on()
		ax3.tick_params(which="both",direction="in")
		# Inset zoom of finer δϕ grid 
		ax4 = fig.add_axes([0.68,0.4,0.2,0.4])
	  ax4.plot(xgrid,xprob,color="purple")
	  ax4.plot(xbin_square,hist_square./maximum(hist_square),color="purple",linewidth=2,alpha=0.5)
	  xlim(2.24,2.37)
	  ax4.minorticks_on()
	  ax4.tick_params(which="both",direction="in",left="false",labelleft="false",
	  		right="true",labelright="true")
	  title=string("IMAGES/likelihoods/",sim,mcmodel,"Moon-",sigma,"secs",nyear,"yrs.png")
	end
  tight_layout()
  savefig(title)
	# clf()
end
