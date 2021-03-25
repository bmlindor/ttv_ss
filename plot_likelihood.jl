using PyPlot,Statistics,Distributions,Optim
rc("font",family="sans-serif")
include("histogram.jl")
function gaussian(x,mu,sig)
  return exp.(-((x .- mu).^2) ./ (2 * sig^.2))
end

function plot_likelihood(jldfit,mcmc,include_moon::Bool=false)
	wide = jldopen("FITS/p3_widefit30.0s40.0yrs.jld2","r")
	p3_wide = 10 .^ range(log10(wide["p3in"]),stop=log10(wide["p3out"]),length=wide["np3"])
	lprob_p3_wide = wide["lprob_p3"]
	pbest_global_wide = wide["pbest_global"]

	xgrid = (10 .^ range(log10(jldfit["p3in"]),stop=log10(jldfit["p3out"]),length=jldfit["np3"])) /365.25
	xprob = exp.((jldfit["lprob_p3"] .-maximum(jldfit["lprob_p3"])))
	bfvalue = jldfit["pbest_global"][12] / 365.25
	param = vec(mcmc["par_mcmc"][:,mcmc["iburn"]:end,12])/365.25
	pbin,phist,pbin_square,hist_square=histogram(param,50)
	label="Period Search Grid [years]"
	color="firebrick"
	if include_moon
		xgrid = range(jldfit["dpin"],stop=jldfit["dpout"],length=jldfit["ndp"])
		xprob = exp.((jldfit["lprob_dp"] .-maximum(jldfit["lprob_dp"])))
		bfvalue = jldfit["pbest_global"][18]
			param = vec(mcmc["par_mcmc"][:,mcmc["iburn"]:end,18])
		pbin,phist,pbin_square,hist_square=histogram(param,50)
		label="Phase Offset Search Grid [radians]"
		color="purple"
	end
	fig = plt.figure(figsize=(8,6))
  ax1 = gca()
	# ax1.axvline(pbest_global_wide[12] /365.25,linestyle="--",color="black",label="Best Fit Value")
	ax1.axvline(11.862615, linestyle=":",color="black",label="True Value ")
	ax1.plot(p3_wide/365.25,exp.((lprob_p3_wide .-maximum(lprob_p3_wide))),color=color) 
	xlabel(label)
	ylabel("Probability")
	ax1.minorticks_on()
	ax1.tick_params(which="both",direction="in")
	ax1.legend(loc="upper left")
	# Inset zoom of finer grid within 
  ax2 = fig.add_axes([0.2,0.3,0.45,0.45])
  # ax2.axvline(bfvalue,linestyle="--",color="black")
  # ax2.errorbar(ttsim1,ttv1,sigtt1,fmt=".",color="black",mec="black",mfc="white")
  ax2.plot(xgrid,xprob,color=color,alpha=0.5)
  ax2.plot(pbin_square,hist_square./maximum(hist_square),color=color,linewidth=2)
  # ax2.plot(xgrid,fit)
  # ax2.set_xlim(4150/365.25,4400/365.25)
  # ax2.minorticks_on()
  ax2.tick_params(which="both",direction="in")
	show()
    # savefig("IMAGES/p3likelihood.png")

# 	if include_moon
# 		clf()
		
# 		lprob_dp = 
# 		pbest_global = jldfit["pbest_global"]
# 		fig = plt.figure(figsize=(7,5))
#         ax3 = gca()
# 		axvline(pbest_global[18],linestyle="--",color="black",label="Best Fit Value")
# 		plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))),color="purple")
# 		# title("Likelihood of Lunar Phase Offset per Earth year")
# 		xlabel("Phase Offset Search Grid [radians]")
# 		ylabel("Likelihood")
# 		ax3.minorticks_on()
# 		ax3.tick_params(which="major",direction="in",top="true",right="true",length=6)
# 		ax3.tick_params(which="minor",direction="in",top="true",right="true",length=2)
# 		ax3.legend(loc="upper left")
#     ax4 = fig.add_axes([0.6,0.6,0.25,0.25])
# #         ax4.axvline(param[18],linestyle="--",color="black")
#     ax4.plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))),color="purple")
#     ax4.minorticks_on()
#     ax4.tick_params(which="minor",direction="in",top="true",right="true",length=2)
#     ax4.tick_params(which="major",direction="in",length=4,
#     left="false",right="false",top="false",bottom="true",
#     labelbottom="true",labeltop="false",labelleft="false",labelright="false")
# 	end
# 	show()
    # savefig("IMAGES/deltaphi.png")
end