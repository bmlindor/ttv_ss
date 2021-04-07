using PyPlot,Statistics,Distributions,Optim
rc("font",family="sans-serif")
include("histogram.jl")
function gaussian(x,mu,sig)
  return exp.(-((x .- mu).^2) ./ (2 * sig^.2))
end

function plot_likelihood(jldfit,mcmc,include_moon::Bool=false)
	# wide = jldopen("FITS/moon_widefit30.0s40.0yrs.jld2","r") # 
	wide = jldopen("FITS/p3_widefit30.0s40.0yrs.jld2","r")
	grid_wide = (10 .^ range(log10(wide["p3in"]),stop=log10(wide["p3out"]),length=wide["np3"])) /365.25
	lprob_wide = exp.((wide["lprob_p3"] .-maximum(wide["lprob_p3"])))
	pbest_global_wide = wide["pbest_global"]

	xgrid = (10 .^ range(log10(jldfit["p3in"]),stop=log10(jldfit["p3out"]),length=jldfit["np3"])) /365.25
	xprob = exp.((jldfit["lprob_p3"] .-maximum(jldfit["lprob_p3"])))
	truex = 11.862615
	lim = 11.4,12.1
	bfvalue = jldfit["pbest_global"][12] /365.25
	param = vec(mcmc["par_mcmc"][:,mcmc["iburn"]:end,12])/365.25
	xbin,xhist,xbin_square,hist_square=histogram(param,50)
	label="Period Search Grid [years]"
	color="firebrick"

	fig = plt.figure(figsize=(8,6))
  ax1 = gca()
	# ax1.axvline(pbest_global_wide[12] /365.25,linestyle="--",color="black",label="Best Fit Value")
	ax1.axvline(truex, linestyle="--",color="black",label=string(truex))
	ax1.plot(grid_wide,lprob_wide,color=color) 
	# ax1.plot(xgrid,xprob,color=color)
	# ax1.legend()
	xlabel(label)
	ylabel("Probability")
	ax1.minorticks_on()
	ax1.tick_params(which="both",direction="in")
	# Inset zoom of finer grid within 
  ax2 = fig.add_axes([0.2,0.4,0.4,0.4])
  # ax2.axvline(bfvalue,linestyle="--",color="black")
  ax2.plot(xgrid,xprob,color=color)
  ax2.plot(xbin_square,hist_square./maximum(hist_square),color=color,linewidth=2,alpha=0.5)
  xlim(lim)
  ax2.minorticks_on()
  ax2.tick_params(which="both",direction="in",left="false",labelleft="false")
	show()
    # savefig("IMAGES/p3likelihood.png")
if include_moon
	# grid_wide = range(wide["dpin"],stop=wide["dpout"],length=wide["ndp"])
	# lprob_wide = exp.((wide["lprob_dp"] .-maximum(wide["lprob_dp"])))
		clf()
	
	xgrid = range(jldfit["dpin"],stop=jldfit["dpout"],length=jldfit["ndp"])
	xprob = exp.((jldfit["lprob_dp"] .-maximum(jldfit["lprob_dp"])))
	truex = 2.31221
	lim = 2.15,2.45
	bfvalue = jldfit["pbest_global"][18]
	param = vec(mcmc["par_mcmc"][:,mcmc["iburn"]:end,18])
	xbin,xhist,xbin_square,hist_square=histogram(param,50)
	label="Phase Offset Search Grid [radians]"
	color="purple"

  fig = plt.figure(figsize=(7,5))
  ax3 = fig.add_subplot(121)
	# ax1.axvline(pbest_global_wide[12] /365.25,linestyle="--",color="black",label="Best Fit Value")
	ax3.axvline(truex, linestyle="--",color="black",label="True Value ")
	ax3.plot(xgrid,xprob,color=color) 
	ax3.plot(xbin_square,hist_square./maximum(hist_square),color=color,linewidth=2,alpha=0.5)
	xlabel(label)
	xlim(-pi/8,3pi)
	ylabel("Probability")
	ax3.minorticks_on()
	ax3.tick_params(which="both",direction="in")
	# ax3.legend(loc="upper left")
	ax4=fig.add_subplot(222)
	# ax4 = fig.add_axes([0.48,0.4,0.4,0.4])
  # ax2.axvline(bfvalue,linestyle="--",color="black")
  ax4.plot(xgrid,xprob,color=color)
  ax4.plot(xbin_square,hist_square./maximum(hist_square),color=color,linewidth=2,alpha=0.5)
  xlim(2.15,2.45)
  ax4.minorticks_on()
  ax4.tick_params(which="both",direction="in",left="false",labelleft="false")
end
# 		ax3.legend(loc="upper left")
#     ax4 = fig.add_axes([0.6,0.6,0.25,0.25])
# #         ax4.axvline(param[18],linestyle="--",color="black")
#     ax4.plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))),color="purple")
end