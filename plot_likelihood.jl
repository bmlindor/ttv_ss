using PyPlot
rc("font",family="serif")



function plot_likelihood(include_moon::Bool=false)
	p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
	deltaphi = range(dpin,stop=dpout,length=ndp)
	nparam=length(pbest_global)
	fig = plt.figure(figsize=(7, 5))
    ax1 = gca()
	ax1.axvline(pbest_global[12]/365.25,linestyle="--",color="black",label="Best Fit Value")
	ax1.plot(p3/365.25,exp.((lprob_p3 .-maximum(lprob_p3))),linewidth=1.25,color="firebrick") 
	# plot!( -5:8,(-5:8).^2,inset = (1,bbox(0.1,0.0,0.4,0.4)),subplot = 2)
	title("Likelihood of Jovian Orbital Period")
	xlabel("Period Search Grid [years]")
	ylabel("Likelihood")
	ax1.minorticks_on()
	ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
	ax1.tick_params(which="minor",direction="in",top="true",right="true")
	ax1.legend(loc="upper left")
    ax2 = fig.add_axes([0.6,0.6,0.25,0.25])
    ax2.axvline(param[12]/365.25,linestyle="--",color="black")
    ax2.plot(p3/365.25,exp.((lprob_p3 .-maximum(lprob_p3))),linewidth=1.25,color="firebrick")
    ax2.set_xlim(4260/365.25,4425/365.25)
    ax2.minorticks_on()
    ax2.tick_params(which="minor",direction="in",top="true",right="true",length=2)
    ax2.tick_params(which="major",direction="in",length=4,
    left="false",right="false",top="false",bottom="true",
    labelbottom="true",labeltop="false",labelleft="false",labelright="false")
	show()
    # savefig("IMAGES/p3likelihood.png")

	if include_moon
		clf()
		fig = plt.figure(figsize=(7, 5))
        ax3 = gca()
		axvline(pbest_global[18],linestyle="--",color="black",label="Best Fit Value")
		plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))),linewidth=1.25,color="purple")
		title("Likelihood of Lunar Phase Offset per Earth year")
		xlabel("Phase Offset Search Grid [radians]")
		ylabel("Likelihood")
		ax3.minorticks_on()
		ax3.tick_params(which="major",direction="in",top="true",right="true",length=6)
		ax3.tick_params(which="minor",direction="in",top="true",right="true",length=2)
		ax3.legend(loc="upper left")
        ax4 = fig.add_axes([0.6,0.6,0.25,0.25])
#         ax4.axvline(param[18]/365.25,linestyle="--",color="black")
        ax4.plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))),linewidth=1.25,color="purple")
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