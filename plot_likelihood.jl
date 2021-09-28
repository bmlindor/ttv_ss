using PyPlot,Statistics,Distributions,Optim,JLD2
rc("font",family="sans-serif")
include("histogram.jl")
avg(x,y) = (x + y)/2
gaussian(x,mu,sig)=exp.(-((x .- mu).^2) ./ (2 * sig^.2))
# Create a plot of likelihood/probability of values in grid
function plot_likelihood(sigma,nyear,sim,model,nbins,include_moon::Bool=false)
	if string(sim)=="EMB" && isfile(string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jl")) 
    mcfile = string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile = string("FITS/fromEMB/",model,"_fit",sigma,"s",nyear,"yrs.jld2")
   else #if isfile(string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")) 
    mcfile = string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile = string("FITS/",model,"_fit",sigma,"s",nyear,"yrs.jld2")
  # else 
  #   return  println("MCMC or FITS file for ",sim," with ",model," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
	wide = jldopen("FITS/wide2_fit30.0s30.0yrs.jld2","r")
	grid_wide = wide["p3"]./365.25
	lprob_wide = exp.((wide["lprob_p3"] .-maximum(wide["lprob_p3"])))
	# EMBfile = string("FITS/fromEMB/p3_fit",sigma,"s",nyear,"yrs.jld2")
	# EMBfit = jldopen(String(EMBfile),"r")
	# mcfile = string("MCMC/fromEMB/p3_mcmc",sigma,"s",nyear,"yrs.jld2")
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
	label="Period Search Grid [years]"
	lim = f["p3in"]/365.25,f["p3out"]/365.25
	color="firebrick"
	fig = plt.figure()
  ax1 = gca()
	# ax1.axvline(pbest_global_wide[ncol] /365.25,linestyle="--",color="black",label="Best Fit Value")
	ax1.axvline(truex,linestyle="--",color="black",label=string(truex))
	ax1.plot(grid,prob,color=color) 
	# ax1.plot(xgrid,xprob,color=color)
	# xlim(minimum(grid_wide),14)
	ax1.legend()
	xlabel(label)
	ylabel("Probability")
	ax1.minorticks_on()
	ax1.tick_params(which="both",direction="in")
	# Inset zoom of finer period grid 
  ax2 = fig.add_axes([0.2,0.4,0.4,0.4])
  # ax2.axvline(bfvalue,linestyle="--",color="black")
  ax2.plot(xgrid,xprob,color=color)
  ax2.plot(xbin_square,hist_square./maximum(hist_square),color=color,linewidth=2,alpha=0.5)
  xlim(lim)
  ax2.minorticks_on()
  ax2.tick_params(which="both",direction="in",left="true",labelleft="true")
	if include_moon
		xgrid = range(f["dpin"],stop=f["dpout"],length=f["ndp"])
		xprob = exp.((f["lprob_dp"] .-maximum(f["lprob_dp"])))
		grid_wide = range(wide["dpin"],stop=wide["dpout"],length=wide["ndp"])
		lprob_wide = exp.((wide["lprob_dp"] .-maximum(wide["lprob_dp"])))
		grid = [grid_wide[1:155];xgrid;grid_wide[198:end]]
		prob = [lprob_wide[1:155];xprob;lprob_wide[198:end]]
		truex = 2.31221
		ncol = 18
		bfvalue = f["best_dp"][ncol]
		param = vec(m["par_mcmc"][:,m["iburn"]:end,ncol])
		xbin,xhist,xbin_square,hist_square=histogram(param,nbins)
		label="Phase Offset Search Grid [radians]"
		lim = f["dpin"],f["dpout"]
		color="purple"
		fig = plt.figure()
	  ax3 = gca()
		ax3.axvline(truex,linestyle="--",color="black",label=string(truex))
		# ax3.plot(grid,prob,color=color) 
		ax3.plot(grid_wide,lprob_wide,color=color) 
		ax3.legend()
		xlabel(label)
		ylabel("Probability")
		ax3.minorticks_on()
		ax3.tick_params(which="both",direction="in")
		# Inset zoom of finer δϕ grid 
		ax4 = fig.add_axes([0.65,0.4,0.2,0.4])
	  ax4.plot(xgrid,xprob,color=color)
	  ax4.plot(xbin_square,hist_square./maximum(hist_square),color=color,linewidth=2,alpha=0.5)
	  xlim(2.23,2.4)
	  ax4.minorticks_on()
	  ax4.tick_params(which="both",direction="in",left="false",labelleft="false",
	  		right="true",labelright="true")
	end
end