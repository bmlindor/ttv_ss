using PyPlot
rc("font",family="sans-serif")

function plot_likelihood(jldfit,include_moon::Bool=false)
	p3 = 10 .^ range(log10(jldfit["p3in"]),stop=log10(jldfit["p3out"]),length=jldfit["np3"])
	lprob_p3 = jldfit["lprob_p3"]
	pbest_global = jldfit["pbest_global"]
	fig = plt.figure(figsize=(8, 6))
  ax1 = gca()
	ax1.axvline(pbest_global[12] /365.25,linestyle="--",color="black",label="Best Fit Value")
	ax1.axvline(11.862615, label="True Value")
	ax1.plot(p3/365.25,exp.((lprob_p3 .-maximum(lprob_p3))),color="firebrick",marker=".") 
	# title("Likelihood of Jovian Orbital Period")
	xlabel("Period Search Grid [years]")
	ylabel("Likelihood")
	ax1.minorticks_on()
	ax1.tick_params(which="both",direction="in")
	# ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
	# ax1.tick_params(which="minor",direction="in",top="true",right="true")
	ax1.legend(loc="upper left")
  ax2 = fig.add_axes([0.2,0.4,0.45,0.35])
  ax2.axvline(pbest_global[12] /365.25,linestyle="--",color="black")
  ax2.axvline(11.862615, label="True Value")
  ax2.plot(p3/365.25,exp.((lprob_p3 .-maximum(lprob_p3))),color="firebrick",marker=".")
  ax2.set_xlim(4150/365.25,4400/365.25)
  ax2.minorticks_on()
  ax2.tick_params(which="both",direction="in",top="true")
  # ax2.tick_params(which="minor",direction="in",top="true",right="true",length=2)
  # ax2.tick_params(which="major",direction="in",length=4,
  # left="false",right="false",top="false",bottom="true",
  # labelbottom="true",labeltop="false",labelleft="false",labelright="false")
	show()
    # savefig("IMAGES/p3likelihood.png")

	if include_moon
		clf()
		deltaphi = range(jldfit["dpin"],stop=jldfit["dpout"],length=jldfit["ndp"])
		lprob_dp = jldfit["lprob_dp"]
		pbest_global = jldfit["pbest_global"]
		fig = plt.figure(figsize=(7, 5))
        ax3 = gca()
		axvline(pbest_global[18],linestyle="--",color="black",label="Best Fit Value")
		plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))),color="purple")
		# title("Likelihood of Lunar Phase Offset per Earth year")
		xlabel("Phase Offset Search Grid [radians]")
		ylabel("Likelihood")
		ax3.minorticks_on()
		ax3.tick_params(which="major",direction="in",top="true",right="true",length=6)
		ax3.tick_params(which="minor",direction="in",top="true",right="true",length=2)
		ax3.legend(loc="upper left")
    ax4 = fig.add_axes([0.6,0.6,0.25,0.25])
#         ax4.axvline(param[18]/365.25,linestyle="--",color="black")
    ax4.plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))),color="purple")
#         ax4.set_xlim(4260/365.25,4425/365.25)
    ax4.minorticks_on()
    ax4.tick_params(which="minor",direction="in",top="true",right="true",length=2)
    ax4.tick_params(which="major",direction="in",length=4,
    left="false",right="false",top="false",bottom="true",
    labelbottom="true",labeltop="false",labelleft="false",labelright="false")
	end
	show()
    # savefig("IMAGES/deltaphi.png")
end