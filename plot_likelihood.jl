using PyPlot,Statistics,Distributions,Optim
rc("font",family="sans-serif")
include("histogram.jl")
avg(x,y) = (x + y)/2
gaussian(x,mu,sig)=exp.(-((x .- mu).^2) ./ (2 * sig^.2))
# Create a plot of likelihood/probability of values in grid
function plot_likelihood(jldfit,mcmc,nbins,include_moon::Bool=false)
	tt,tt0,sigtt,ttmodel = jldfit["tt"],jldfit["tt0"],jldfit["sigtt"],jldfit["ttmodel"]
  pbest_global = jldfit["pbest_global"]
  nplanet,ntrans = jldfit["nplanet"],jldfit["ntrans"]
	par_mcmc = mcmc["par_mcmc"]
  lprob_mcmc = mcmc["lprob_mcmc"]
  param = mcmc["param"]
  iburn = mcmc["iburn"]
  nwalkers = mcmc["nwalkers"]
  nsteps = mcmc["nsteps"]

	wide = jldopen("FITS/p3_widefit30.0s30.0yrs.jld2","r")
	grid_wide = (10 .^ range(log10(wide["p3in"]),stop=log10(wide["p3out"]),length=wide["np3"])) /365.25
	lprob_wide = exp.((wide["lprob_p3"] .-maximum(jldfit["lprob_p3"])))
	pbest_global_wide = wide["pbest_global"]
	xgrid = (10 .^ range(log10(jldfit["p3in"]),stop=log10(jldfit["p3out"]),length=jldfit["np3"])) /365.25
	xprob = exp.((jldfit["lprob_p3"] .-maximum(jldfit["lprob_p3"])))
	truex = 11.862615
	bfvalue = jldfit["pbest_global"][12] /365.25
	param = vec(mcmc["par_mcmc"][:,mcmc["iburn"]:end,12])/365.25
	xbin,xhist,xbin_square,hist_square=histogram(param,nbins)
	label="Period Search Grid [years]"
	lim = jldfit["p3in"]/365.25,jldfit["p3out"]/365.25
	color="firebrick"
	println("Simulated with σ= ",sigtt[1]*3600*24," second noise")

	fig = plt.figure()
  ax1 = gca()
	# ax1.axvline(pbest_global_wide[12] /365.25,linestyle="--",color="black",label="Best Fit Value")
	ax1.axvline(truex,linestyle="--",color="black",label=string(truex))
	grid = [grid_wide[1:end-22];avg(grid_wide[end-21],xgrid[1]);xgrid]#;grid_wide[end-1:end]]
	prob = [lprob_wide[1:end-22];avg(lprob_wide[end-21],xprob[1]);xprob]#;lprob_wide[end-1:end]]
	ax1.plot(grid,prob,color=color) 
	# ax1.plot(xgrid,xprob,color=color)
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
	show()
    # savefig("IMAGES/p3likelihood.png")
if include_moon
		clf()
	wide = jldopen("FITS/moon_widefit30.0s40.0yrs.jld2","r") 
	grid_wide = range(wide["dpin"],stop=wide["dpout"],length=wide["ndp"])
	lprob_wide = exp.((wide["lprob_dp"] .-maximum(wide["lprob_dp"])))
	xgrid = range(jldfit["dpin"],stop=jldfit["dpout"],length=jldfit["ndp"])
	xprob = exp.((jldfit["lprob_dp"] .-maximum(jldfit["lprob_dp"])))
	truex = 2.31221
	bfvalue = jldfit["pbest_global"][18]
	param = vec(mcmc["par_mcmc"][:,mcmc["iburn"]:end,18])
	xbin,xhist,xbin_square,hist_square=histogram(param,nbins)
	label="Phase Offset Search Grid [radians]"
	lim = jldfit["dpin"],jldfit["dpout"]
	color="purple"
	# println(grid_wide[end-31:end])
  fig = plt.figure()
  ax3 = gca()
	ax3.axvline(truex,linestyle="--",color="black",label=string(truex))
	grid = [grid_wide[1:end-50];avg(grid_wide[end-49],xgrid[1]);xgrid;grid_wide[end-30:end]]
	prob = [lprob_wide[1:end-50];avg(lprob_wide[end-49],xprob[1]);xprob;lprob_wide[end-30:end]]
	ax3.plot(grid,prob,color=color) 
	# ax3.plot(grid_wide,lprob_wide,color=color) 
	# ax3.plot(xbin_square,hist_square./maximum(hist_square),color=color,linewidth=2,alpha=0.5)
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