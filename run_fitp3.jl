# Julia v1.1
using JLD2, DelimitedFiles, PyPlot
include("sim_times.jl")
include("fit_mysteryplanet3.jl")
# sig_grid = [10.0, 15.0, 30.0, 45.0, 60.0]
# p3in = 500.0; p3out = 18000.0; np3 = 1000
p3in = 4230.0; p3out = 4400.0; np3 = 10
p3 = 10 .^ range(log10(p3in),stop=log10(p3out),length=np3)
function run_fitp3(p3in, p3out, np3, sigma)
	# sim = sim_times(2.4332825e6, 2.4515445e6, 1000, true, sigma)
	file = string("INPUTS/tt_data",sigma,"s.txt")
	fit = fit_mysteryplanet3(file,p3in,p3out,np3,6,true, sigma)
end

run_fitp3(p3in, p3out, np3, 30.0)
# function make_plot(filename::String)
# 	# chi_best, pbest, chi_p3, param_p3
# 	data = readdlm(filename)
# 	chi_best, pbest, chi_p3, param_p3 = data[1], data[2], data[3], data[4]
# 	plot(p3/365.25,exp.(-0.5*(chi_p3 .-minimum(chi_p3)))) #to show that max likelihood peaks at actual period
#     xlabel("Period of planet 3 [years]")
#     ylabel("Likelihood")
# end