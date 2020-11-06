using PyPlot
rc("font",family="serif")

function plot_likelihood(include_moon::Bool=false)
p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
deltaphi = range(dpin,stop=dpout,length=ndp)
nparam=length(pbest_global)

figsize=(8,6)
if include_moon
	plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))),color="purple")
	xlabel(L"$\Delta \phi$ Grid [radians]")
	name = string("IMAGES/deltaphi.png")

else
	plot(p3/365.25,exp.((lprob_p3 .-maximum(lprob_p3))),color="firebrick") 
	# plot!( -5:8,(-5:8).^2,inset = (1,bbox(0.1,0.0,0.4,0.4)),subplot = 2)
	xlabel("Period Grid [years]")
	name = string("IMAGES/p3likelihood.png")
end
ylabel("Likelihood")
tick_params(direction="in")
# savefig(name)
end