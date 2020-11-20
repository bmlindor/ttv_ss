using PyPlot
rc("font",family="serif")



function plot_likelihood(include_moon::Bool=false)
	p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
	deltaphi = range(dpin,stop=dpout,length=ndp)
	nparam=length(pbest_global)
	fig, ax1 = subplots(figsize=(6,4))
	axvline(pbest_global[12]/365.25,linestyle="--",color="black",label="Best Fit Value")
	plot(p3/365.25,exp.((lprob_p3 .-maximum(lprob_p3))),linewidth=1.25,color="firebrick") 
	# plot!( -5:8,(-5:8).^2,inset = (1,bbox(0.1,0.0,0.4,0.4)),subplot = 2)
	title("Likelihood of Jovian Orbital Period")
	xlabel("Period Search Grid [years]")
	ylabel("Likelihood")
	ax1.minorticks_on()
	ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
	ax1.tick_params(which="minor",direction="in",top="true",right="true")
	legend()
	tight_layout()
	show()
    # savefig("IMAGES/p3likelihood.png")

	if include_moon
		clf()
		fig, ax2 = subplots(figsize=(6,4))
		axvline(pbest_global[18],linestyle="--",color="black",label="Best Fit Value")
		plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))),linewidth=1.25,color="purple")
		title("Likelihood of Lunar Phase Offset per Earth year")
		xlabel("Phase Offset Search Grid [radians]")
		ylabel("Likelihood")
		ax2.minorticks_on()
		ax2.tick_params(which="major",direction="in",top="true",right="true",length=6)
		ax2.tick_params(which="minor",direction="in",top="true",right="true",length=2)
		legend()
		tight_layout()
	end
	show()
    # savefig("IMAGES/deltaphi.png")
end