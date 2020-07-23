# Julia v1.1
using JLD2, DelimitedFiles, PyPlot
include("sim_times.jl")
include("fit_mysteryplanet3.jl")
function run_fitp3()
	sim = sim_times(2.4332825e6, 2.4515445e6, 1000, true, 30.0)
	fit = fit_mysteryplanet3("noisy_ttvenus.txt","noisy_ttearth.txt",4230.0,4400.0,100,6,true)
end
# function make_plot(filename::String)
# 	# chi_best, pbest, chi_p3, param_p3
# 	data = readdlm(filename)
# 	chi_best, pbest, chi_p3, param_p3 = data[1], data[2], data[3], data[4]
# 	plot(p3/365.25,exp.(-0.5*(chi_p3 .-minimum(chi_p3)))) #to show that max likelihood peaks at actual period
#     xlabel("Period of planet 3 [years]")
#     ylabel("Likelihood")
# end