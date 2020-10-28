using PyPlot
rc("font",family="serif")

function plot_likelihood(include_moon::Bool=false)

figsize=(8,6)
if include_moon
	plot(deltaphi,exp.((lprob_dp .-maximum(lprob_dp))))
	xlabel("δϕ of Moon [radians]")
	name = string("IMAGES/deltaphi.png")

else
	plot(p3/365.25,exp.((lprob_p3 .-maximum(lprob_p3)))) 
	# plot!( -5:8,(-5:8).^2,inset = (1,bbox(0.1,0.0,0.4,0.4)),subplot = 2)
	xlabel("Period of planet 3 [years]")
	name = string("IMAGES/p3likelihood.png")
end
ylabel("Likelihood")
# savefig(name)
end