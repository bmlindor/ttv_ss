using PyPlot,Statistics,Distributions,JLD2
rc("font",family="sans-serif")
include("histogram.jl")

function plot_profile(xgrid::Array{Float64,1},lprob::Array{Float64,1},truex::Float64,nbins::Int,color::String,pname::String,sim_obs_label::String,label_xloc::Real)
	fig=figure(figsize=(6,6))
 	subplots_adjust(hspace=0.05,wspace=0.05)
	  #ax1=gca()
  xprob=exp.(lprob .- maximum(lprob))
	xgrid_in_yrs =xgrid ./365.25 
	#scatter(xgrid_in_yrs,xprob)
	plot(xgrid_in_yrs,xprob,color=color,label=sim_obs_label) 
	#xbin,xhist,xbin_square,hist_square=histogram(xgrid_in_yrs,nbins)
	#plot(xbin_square,hist_square./maximum(hist_square),linewidth=2,alpha=0.5)
	axvline(truex,linestyle="--",color="black")
	text(truex + truex/100,1.01,pname)
	text(label_xloc,1.05,sim_obs_label)
	#legend(loc="upper left")
	xlabel("Planet Period Search Grid [years]")
	xlim(minimum(xgrid_in_yrs)-minimum(xgrid_in_yrs)/10,maximum(xgrid_in_yrs)+.05)
	ylabel("Probability")
	ylim(0,1.2)
	minorticks_on()
	tick_params(which="both",direction="in")
	#tight_layout()
end
# Create likelihood/probability plots of search grids
function plot_grid(sigma,nyear,grid_type_nplanet,per_col,true_per,color,pname,case_num,label_xloc)
	if case_num==1
	file=string("/astro/users/blindor/research/ttv_ss/grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 1"
	elseif case_num==2
	file=string("/astro/users/blindor/research/ttv_ss/grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 2"
	end
	data,header=readdlm(file,',',header=true)
	sim_obs_label= string(case_label," [",nyear," yr span]",'\n',L"$\sigma_{obs}=$",sigma," sec")
 plot_profile(data[:,per_col],data[:,end],true_per,100,color,pname,sim_obs_label,label_xloc)
	savefig("IMAGES/case",case_num,"_",grid_type_nplanet,pname,sigma,"secs",year,"yrs.png")
end

function plot_prob(sigma::Real,nyear::Real,obs::String,fitmodel::String,mcmodel::String,nbins::Int,include_moon::Bool=false)
	if obs=="fromEMB" && isfile(string("MCMC/fromEMB/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")) 
    mcfile=string("MCMC/fromEMB/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile=string("FITS/fromEMB/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
  elseif isfile(string("MCMC/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")) 
    mcfile=string("MCMC/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile=string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
    println("Loading...",mcfile," and ",fitfile)
  else 
    return  println("MCMC or FITS file for ",sim," with ",mcmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
 	m=jldopen(String(mcfile),"r")
	f=jldopen(String(fitfile),"r")

	# p3=vec(f["p3"])./365.25
	# p3prob=exp.((vec(f["lprob_p3"]) .-maximum(vec(f["lprob_p3"]))))
	# p3_cur=11.862615
	# bfvalue=f["best_p3"][12] /365.25	
	# param=vec(m["par_mcmc"][:,m["iburn"]:end,12])./365.25
	# prob(p3,p3prob,p3_cur,nbins,param,"firebrick","Jupiter")
	# title=string("IMAGES/likelihoods/",sim,mcmodel,"planet3-",sigma,"secs",nyear,"yrs.png")
	# show()
	# savefig(title)
 #  if fitmodel=="p4"
		fig=figure(figsize=(6,6))
  	subplots_adjust(hspace=0.05,wspace=0.05)
	  ax1=gca()
	  xgrid4=f["p4"]./365.25
		xprob4=exp.((f["lprob_p4"] .-maximum(f["lprob_p4"])))
		# grid=[grid_wide[1:55];xgrid;grid_wide[60:end]]
		# prob=[lprob_wide[1:55];xprob;lprob_wide[60:end]]
		bfvalue4=f["best_p4"][12] /365.25	
		truex4=1.88
		param4=vec(m["par_mcmc"][:,m["iburn"]:end,12])./365.25
		xbin4,xhist4,xbin_square4,hist_square4=histogram(param4,nbins)
		ax1.axvline(truex4,linestyle="--",color="black",label=string(truex4))
		ax1.plot(xgrid4,xprob4,color="darkcyan")
		ax1.text(2,1,"Mars")
		ax1.set_xlim(1,5)
	#   title=string("IMAGES/likelihoods/",sim,mcmodel,"planet4-",sigma,"secs",nyear,"yrs.png")
	# else
	# 	ax1.plot(xgrid,xprob,color="firebrick") 
	# 	# Inset zoom of finer period grid 
	#  #  ax2=fig.add_axes([0.2,0.4,0.4,0.5])
	#  #  # ax2.axvline(bfvalue,linestyle="--",color="black")
	#  #  ax2.plot(xgrid,xprob,color="firebrick")
	#  #  ax2.plot(xbin_square,hist_square./maximum(hist_square),color="firebrick",linewidth=2,alpha=0.5)
	#  #  xlim(lim)
	# 	# # ax2.minorticks_on()
	#  #  ax2.tick_params(which="both",direction="in",left="true",labelleft="true")
	# 	savefig(title)
	#   tight_layout()
 #  end
		show()
end
# Plot wide likelihood of values in search grids, with probability inset
function plot_likelihood(sigma::Real,nyear::Real,sim::String,fitmodel::String,mcmodel::String,nbins::Int,include_moon::Bool=false)
	if String(sim)=="EMB" && isfile(string("MCMC/fromEMB/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")) 
    mcfile=string("MCMC/fromEMB/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile=string("FITS/fromEMB/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
  elseif isfile(string("MCMC/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")) 
    mcfile=string("MCMC/",mcmodel,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile=string("FITS/",fitmodel,"_fit",sigma,"s",nyear,"yrs.jld2")
    println("Loading...",mcfile," and ",fitfile)
  else 
    return  println("MCMC or FITS file for ",sim," with ",mcmodel," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
 
	wide=jldopen("FITS/wide2_fit30.0s30.0yrs.jld2","r")
	grid_wide=wide["p3"]./365.25
	lprob_wide=exp.((wide["lprob_p3"] .-maximum(wide["lprob_p3"])))
	m=jldopen(String(mcfile),"r")
	f=jldopen(String(fitfile),"r")

	xgrid=f["p3"]./365.25
	xprob=exp.((f["lprob_p3"] .-maximum(f["lprob_p3"])))
	grid=[grid_wide[1:155];xgrid;grid_wide[198:end]]
	prob=[lprob_wide[1:155];xprob;lprob_wide[198:end]]
	# println("Simulated with σ= ",sigma," second noise, for ",nyear," observation years.")
	truex=11.862615
	ncol=12
	bfvalue=f["best_p3"][ncol] /365.25	
	param=vec(m["par_mcmc"][:,m["iburn"]:end,ncol])./365.25
	xbin,xhist,xbin_square,hist_square=histogram(param,nbins)
	label="Planet Period Search Grid [years]"
	lim=minimum(xgrid),maximum(xgrid)
	fig=plt.figure()
  ax1=gca()
	ax1.axvline(truex,linestyle="--",color="black",label=string(truex))
	ax1.text(12,1,"Jupiter")
	# ax1.plot(grid_wide,lprob_wide,linestyle="-",color="grey")
	ax1.plot(grid,prob,color="black",alpha=0.75)
	# xlim(minimum(grid_wide),14)
	# ax1.legend(loc="upper left")
	xlabel(label)
	ylabel("Probability")
	ax1.minorticks_on()
	ax1.tick_params(which="both",direction="in")
  # savefig(title)
  # clf()
  if fitmodel=="p4"
	  xgrid4=f["p4"]./365.25
		xprob4=exp.((f["lprob_p4"] .-maximum(f["lprob_p4"])))
		# grid=[grid_wide[1:55];xgrid;grid_wide[60:end]]
		# prob=[lprob_wide[1:55];xprob;lprob_wide[60:end]]
		bfvalue4=f["best_p4"][ncol] /365.25	
		truex4=1.88
		param4=vec(m["par_mcmc"][:,m["iburn"]:end,ncol])./365.25
		xbin4,xhist4,xbin_square4,hist_square4=histogram(param4,nbins)
		ax1.plot(xgrid,xprob,color="firebrick") 
		ax1.axvline(truex4,linestyle="--",color="black",label=string(truex4))
		ax1.plot(xgrid4,xprob4,color="darkcyan")
		ax1.text(2,1,"Mars")
		ax1.set_xlim(1,14.3)
	  title=string("IMAGES/likelihoods/",sim,mcmodel,"Mars-",sigma,"secs",nyear,"yrs.png")
	else
		ax1.plot(grid,prob,color="firebrick") 
		# Inset zoom of finer period grid 
	  ax2=fig.add_axes([0.2,0.4,0.4,0.5])
	  # ax2.axvline(bfvalue,linestyle="--",color="black")
	  ax2.plot(xgrid,xprob,color="firebrick")
	  ax2.plot(xbin_square,hist_square./maximum(hist_square),color="firebrick",linewidth=2,alpha=0.5)
	  xlim(lim)
		# ax2.minorticks_on()
	  ax2.tick_params(which="both",direction="in",left="true",labelleft="true")
	  title=string("IMAGES/likelihoods/",sim,mcmodel,"Jupiter-",sigma,"secs",nyear,"yrs.png")
	  tight_layout()
  end
	if include_moon
		xgrid=range(f["dpin"],stop=f["dpout"],length=f["ndp"])
		xprob=exp.((f["lprob_dp"] .-maximum(f["lprob_dp"])))
		grid_wide=range(wide["dpin"],stop=wide["dpout"],length=wide["ndp"])
		lprob_wide=exp.((wide["lprob_dp"] .-maximum(wide["lprob_dp"])))
		grid=[grid_wide[1:60];xgrid;grid_wide[75:end]]
		prob=[lprob_wide[1:60];xprob;lprob_wide[75:end]]
		truex=2.31586
		ncol=18
		bfvalue=f["best_dp"][ncol]
		param=vec(m["par_mcmc"][:,m["iburn"]:end,ncol])
		xbin,xhist,xbin_square,hist_square=histogram(param,nbins)
		label="Satellite Phase Offset Search Grid [radians]"
		lim=f["dpin"],f["dpout"]
		fig=plt.figure()
	  ax3=gca()
		ax3.axvline(truex,linestyle="--",color="black",label=string(truex))
		ax3.text(2.4,1,"Moon")
		ax3.plot(grid,prob,color="purple") 
		# tight_layout()
		# ax3.plot(grid_wide,lprob_wide,color="red") 
		xlabel(label)
		ylabel("Probability")
		ax3.minorticks_on()
		ax3.tick_params(which="both",direction="in")
		# Inset zoom of finer δϕ grid 
		ax4=fig.add_axes([0.68,0.4,0.2,0.4])
	  ax4.plot(xgrid,xprob,color="purple")
	  ax4.plot(xbin_square,hist_square./maximum(hist_square),color="purple",linewidth=2,alpha=0.5)
	  xlim(2.24,2.37)
	  ax4.minorticks_on()
	  ax4.tick_params(which="both",direction="in",left="false",labelleft="false",
	  		right="true",labelright="true")
	  title=string("IMAGES/likelihoods/",sim,mcmodel,"Moon-",sigma,"secs",nyear,"yrs.png")
	end
  tight_layout()
  # savefig(title)
	show()
end
