using PyPlot,Statistics,Distributions,JLD2,DelimitedFiles
rc("font",family="sans-serif")
rc("lines",linewidth=2)
include("histogram.jl")
xprob(lprob)=exp.(lprob .- maximum(lprob))
function plot_profile(xgrid::Array{Float64,1},lprob::Array{Float64,1},truex::Float64,nbins::Int,color::String,pname::String,sim_obs_label::String)

	  #ax1=gca()
  xprob=exp.(lprob .- maximum(lprob))
	xgrid_in_yrs =xgrid ./365.25 
	#scatter(xgrid_in_yrs,xprob)
	if color=="salmon"
	plot(xgrid_in_yrs,xprob,color=color) 
	else
	plot(xgrid_in_yrs,xprob,label=sim_obs_label)
	end
	#xbin,xhist,xbin_square,hist_square=histogram(xgrid_in_yrs,nbins)
	#plot(xbin_square,hist_square./maximum(hist_square),alpha=0.5)
end
# Create likelihood/probability plots of search grids
function per_grid(sigma,nyear,grid_type_nplanet,per_col,true_per,color,pname,case_num,label_xloc)
	if case_num==1
	file=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 1"
	elseif case_num==2
	file=string("grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 2"
	end
	data,header=readdlm(file,',',header=true)
	sim_obs_label= string(case_label," [",nyear," yr span]",'\n',L"$\sigma_{obs}=$",sigma," sec")
	# save_as =string("IMAGES/wide_grids/case",case_num,"_",grid_type_nplanet,pname,sigma,"secs",nyear,"yrs.png")
	axvline(true_per,linestyle="--",color="black")
	text(true_per + true_per/100,1.01,pname,fontsize="large")
	text(label_xloc,1.05,sim_obs_label,fontsize="medium")
 	plot(data[:,per_col]./365.25,xprob(data[:,end]),color)
 	xlabel("Planet Period Search Grid",fontsize="x-large")
 	ylabel("Probability",fontsize="x-large")
 	ylim(0,1.2)
end
function moon_grid(sigma,nyear,grid_type_nplanet,per_col,true_per,color,pname,case_num,label_xloc)
	if case_num==1
	file=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 1"
	elseif case_num==2
	file=string("grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 2"
	end
	data,header=readdlm(file,',',header=true)
	sim_obs_label= string(case_label," [",nyear," yr span]",'\n',L"$\sigma_{obs}=$",sigma," sec")
	# save_as =string("IMAGES/wide_grids/case",case_num,"_",grid_type_nplanet,pname,sigma,"secs",nyear,"yrs.png")
	axvline(true_per,linestyle="--",color="black")
	text(true_per + true_per/100,1.01,pname,fontsize="large")
	text(label_xloc,1.05,sim_obs_label,fontsize="medium")
 	plot(data[:,per_col],xprob(data[:,end]),color)
 	xlabel("Planet Period Search Grid",fontsize="x-large")
 	ylabel("Probability",fontsize="x-large")
 	ylim(0,1.2)
end

function compare_grid(nyear,grid_type_nplanet,per_col,true_per,color,pname,case_num)
	sigma=30
	if case_num==1
	file=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
		case_label=string("Case 1"," [",nyear," yr span]")
	elseif case_num==2
	file=string("grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
		case_label=string("Case 2"," [",nyear," yr span]")
	end
	data,header=readdlm(file,',',header=true)
	sim_obs_label= string(case_label," [",nyear," yr span]",'\n',L"$\sigma_{obs}=$",sigma," sec")
	# save_as =string("IMAGES/wide_grids/case",case_num,"_",grid_type_nplanet,pname,sigma,"secs",nyear,"yrs.png")
	axvline(true_per,linestyle="--",color="black")
	text(true_per + true_per/100,1.01,pname,fontsize="large")
 	plot(data[:,per_col]./365.25,xprob(data[:,end]),color=color,label=case_label)
 	xlabel("Planet Period Search Grid",fontsize="x-large")
 	ylabel("Probability",fontsize="x-large")
 	# legend()
 	ylim(0,1.2)
end
function plot_grid(sigma,nyear,grid_type_nplanet,per_col,true_per,pname,case_num,label_xloc)
	if case_num==1
	file=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 1"

	elseif case_num==2
	file=string("grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 2"
		# case_label=string("Case 1"," [",nyear," yr span]")
	end
	data,header=readdlm(file,',',header=true)
	sim_obs_label= string(case_label," [",nyear," yr span]",'\n',L"$\sigma_{obs}=$",sigma," sec")
	axvline(true_per,linestyle="--",color="black")
	text(true_per + true_per/100,1.01,pname)
	text(label_xloc,1.05,sim_obs_label)
 	plot(data[:,per_col],xprob(data[:,end]))
 	legend()
 	ylim(0,1.2)
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
	#  #  ax2.plot(xbin_square,hist_square./maximum(hist_square),color="firebrick",alpha=0.5)
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
	  ax2.plot(xbin_square,hist_square./maximum(hist_square),color="firebrick",alpha=0.5)
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
	  ax4.plot(xbin_square,hist_square./maximum(hist_square),color="purple",alpha=0.5)
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
function compare_sigs(nyear,grid_type_nplanet,per_col,true_per,pname,case_num)
	if case_num==1
	file1=string("grid/fromEMB/",grid_type_nplanet,"_grid10s",nyear,"yrs.csv")	
	file2=string("grid/fromEMB/",grid_type_nplanet,"_grid30s",nyear,"yrs.csv")
	file3=string("grid/fromEMB/",grid_type_nplanet,"_grid60s",nyear,"yrs.csv")
	file4=string("grid/fromEMB/",grid_type_nplanet,"_grid80s",nyear,"yrs.csv")
	file5=string("grid/fromEMB/",grid_type_nplanet,"_grid90s",nyear,"yrs.csv")
	# file6=string("grid/fromEMB/",grid_type_nplanet,"_grid100s",nyear,"yrs.csv")
	# file7=string("grid/fromEMB/",grid_type_nplanet,"_grid110s",nyear,"yrs.csv")
	# file7=string("grid/fromEMB/",grid_type_nplanet,"_grid120s",nyear,"yrs.csv")
	case_label=string("Case 1"," [",nyear," yr span]")
	elseif case_num==2
	file1=string("grid/",grid_type_nplanet,"_grid10s",nyear,"yrs.csv")
	file2=string("grid/",grid_type_nplanet,"_grid30s",nyear,"yrs.csv")
	file3=string("grid/",grid_type_nplanet,"_grid60s",nyear,"yrs.csv")
	file4=string("grid/",grid_type_nplanet,"_grid80s",nyear,"yrs.csv")
	file4=string("grid/",grid_type_nplanet,"_grid90s",nyear,"yrs.csv")
	case_label=string("Case 2"," [",nyear," yr span]")
	end
	data1,header1=readdlm(file1,',',header=true)
	data2,header2=readdlm(file2,',',header=true)
	data3,header3=readdlm(file3,',',header=true)
	data4,header4=readdlm(file4,',',header=true)
		# save_as =string("IMAGES/wide_grids/case",case_num,"_",grid_type_nplanet,pname,sigma,"secs",nyear,"yrs.png")
fig,ax=subplots(figsize=(5,5))#,dpi=150)
ax.plot(data1[:,per_col]./365.25 ,xprob(data1[:,end]),label=L"$\sigma_{obs}=$ 10 sec")
ax.plot(data2[:,per_col]./365.25 ,xprob(data2[:,end]),label="30 sec")
ax.plot(data3[:,per_col]./365.25 ,xprob(data3[:,end]),label="60 sec",linestyle="--")
ax.plot(data4[:,per_col]./365.25 ,xprob(data4[:,end]),label="80 sec",linestyle="-.")
axvline(true_per,linestyle="--",color="black")
text(true_per + true_per/100,1.01,pname,fontsize="large")
# text(label_xloc,1.05,sim_obs_label)
# legend(loc="upper right",fontsize="medium",title=case_label,title_fontsize="medium",bbox_to_anchor=(0.,1.02,1.,.102),ncol=4,mode="expand",borderaxespad=0.0)
legend(fontsize="medium",title_fontsize="medium",title="case_label")
xlabel("Planet Period Search Grid [years]",fontsize="x-large")
# xlim(minimum(xgrid_in_yrs)-minimum(xgrid_in_yrs)/10,maximum(xgrid_in_yrs)+.05)
ylabel("Probability",fontsize="x-large")
ylim(0,1.2)
minorticks_on()
tick_params(which="both",direction="in")
tight_layout()
end

function compare_yrs(sigma,grid_type_nplanet,per_col,true_per,pname,case_num)
	if case_num==1
	file1=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s15yrs.csv")	
	file2=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s16yrs.csv")
	file3=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s17yrs.csv")
	file4=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s18yrs.csv")
	file5=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s19yrs.csv")
	file6=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s20yrs.csv")
	file7=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s21yrs.csv")
	file8=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s22yrs.csv")
	file9=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s23yrs.csv")
	file10=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s24yrs.csv")
	case_label=string("Case 1"," [σ=",sigma," sec]")
	elseif case_num==2
	file1=string("grid/",grid_type_nplanet,"_grid",sigma,"s15yrs.csv")
	file2=string("grid/",grid_type_nplanet,"_grid",sigma,"s16yrs.csv")
	file3=string("grid/",grid_type_nplanet,"_grid",sigma,"s17yrs.csv")
	file4=string("grid/",grid_type_nplanet,"_grid",sigma,"s18yrs.csv")
	file5=string("grid/",grid_type_nplanet,"_grid",sigma,"s19yrs.csv")
	file6=string("grid/",grid_type_nplanet,"_grid",sigma,"s20yrs.csv")
	file7=string("grid/",grid_type_nplanet,"_grid",sigma,"s21yrs.csv")
	file8=string("grid/",grid_type_nplanet,"_grid",sigma,"s22yrs.csv")
	file9=string("grid/",grid_type_nplanet,"_grid",sigma,"s23yrs.csv")
	file10=string("grid/",grid_type_nplanet,"_grid",sigma,"s24yrs.csv")
	case_label=string("Case 2"," [σ=",sigma," sec]")
	end
	data1,header1=readdlm(file1,',',header=true)
	data2,header2=readdlm(file2,',',header=true)
	data3,header3=readdlm(file3,',',header=true)
	data4,header4=readdlm(file4,',',header=true)
	data5,header5=readdlm(file5,',',header=true)
	data6,header6=readdlm(file6,',',header=true)
	data7,header7=readdlm(file7,',',header=true)
	data8,header8=readdlm(file8,',',header=true)
	data9,header9=readdlm(file9,',',header=true)
	data10,header10=readdlm(file10,',',header=true)
fig,ax=subplots(figsize=(8,6))#,dpi=150)
title("Search for 3rd Planet (i.e. Jupiter)",fontsize="xx-large")
ax.plot(data1[:,per_col]./365.25 ,xprob(data1[:,end]),label="15 yrs")
# ax.plot(data2[:,per_col]./365.25 ,xprob(data2[:,end]),label="16 yrs")
ax.plot(data3[:,per_col]./365.25 ,xprob(data3[:,end]),label="17 yrs")
# ax.plot(data4[:,per_col]./365.25 ,xprob(data4[:,end]),label="18 yrs")
ax.plot(data5[:,per_col]./365.25 ,xprob(data5[:,end]),label="19 yrs")
# ax.plot(data6[:,per_col]./365.25 ,xprob(data6[:,end]),label="20 yrs")
ax.plot(data7[:,per_col]./365.25 ,xprob(data7[:,end]),label="21 yrs")
# ax.plot(data8[:,per_col]./365.25 ,xprob(data8[:,end]),label="22 yrs")
ax.plot(data9[:,per_col]./365.25 ,xprob(data9[:,end]),label="23 yrs")
# ax.plot(data10[:,per_col]./365.25 ,xprob(data10[:,end]),label="24 yrs")
axvline(true_per,linestyle="--",color="black")
text(true_per + true_per/100,1.01,pname,fontsize="large")
#inset 
# text(label_xloc,1.05,sim_obs_label)
legend(loc="upper right",fontsize="large",title=case_label,title_fontsize="large",bbox_to_anchor=(0.,0.9,1.,.102),ncol=5,mode="expand",borderaxespad=0.0)
# legend(fontsize="medium",title_fontsize="medium",title=case_label)
if grid_type_nplanet=="p4"
ax2=fig.add_axes([0.5,0.2,0.3,0.5])
ax2.plot(data1[:,per_col]./365.25 ,xprob(data1[:,end]),label="15 yrs")
# ax2.plot(data2[:,per_col]./365.25 ,xprob(data2[:,end]),label="16 yrs")
ax2.plot(data3[:,per_col]./365.25 ,xprob(data3[:,end]),label="17 yrs")
# ax2.plot(data4[:,per_col]./365.25 ,xprob(data4[:,end]),label="18 yrs")
ax2.plot(data5[:,per_col]./365.25 ,xprob(data5[:,end]),label="19 yrs")
# ax2.plot(data6[:,per_col]./365.25 ,xprob(data6[:,end]),label="20 yrs")
ax2.plot(data7[:,per_col]./365.25 ,xprob(data7[:,end]),label="21 yrs")
# ax2.plot(data8[:,per_col]./365.25 ,xprob(data8[:,end]),label="22 yrs")
ax2.plot(data9[:,per_col]./365.25 ,xprob(data9[:,end]),label="23 yrs")
# ax2.plot(data10[:,per_col]./365.25 ,xprob(data10[:,end]),label="24 yrs")
ax2.set_xlim(1.81,1.9)
ax2.axvline(true_per,linestyle="--",color="black")
ax2.tick_params(which="both",direction="in")
end
ax.set_xlabel("Planet Period Search Grid [years]",fontsize="xx-large")
# xlim(minimum(xgrid_in_yrs)-minimum(xgrid_in_yrs)/10,maximum(xgrid_in_yrs)+.05)
ax.set_ylabel("Probability",fontsize="xx-large")
ax.set_ylim(0,1.2)
ax.minorticks_on()
tight_layout()
end
function wide_yrs(sigma,grid_type_nplanet,per_col,true_per,pname,case_num)
	if case_num==1
	file2=string("grid_wide/fromEMB/",grid_type_nplanet,"_grid",sigma,"s16yrs.csv")
	file4=string("grid_wide/fromEMB/",grid_type_nplanet,"_grid",sigma,"s18yrs.csv")
	file6=string("grid_wide/fromEMB/",grid_type_nplanet,"_grid",sigma,"s20yrs.csv")
	file8=string("grid_wide/fromEMB/",grid_type_nplanet,"_grid",sigma,"s22yrs.csv")
	file10=string("grid_wide/fromEMB/",grid_type_nplanet,"_grid",sigma,"s24yrs.csv")
	case_label=string("Case 1"," [σ=",sigma," sec]")
	elseif case_num==2
	file2=string("grid_wide/",grid_type_nplanet,"_grid",sigma,"s16yrs.csv")
	file4=string("grid_wide/",grid_type_nplanet,"_grid",sigma,"s18yrs.csv")
	file6=string("grid_wide/",grid_type_nplanet,"_grid",sigma,"s20yrs.csv")
	file8=string("grid_wide/",grid_type_nplanet,"_grid",sigma,"s22yrs.csv")
	file10=string("grid_wide/",grid_type_nplanet,"_grid",sigma,"s24yrs.csv")
	case_label=string("Case 2"," [σ=",sigma," sec]")
	end
	data2,header2=readdlm(file2,',',header=true)
	data4,header4=readdlm(file4,',',header=true)
	data6,header6=readdlm(file6,',',header=true)
	data8,header8=readdlm(file8,',',header=true)
	data10,header10=readdlm(file10,',',header=true)
	fig,ax=subplots(figsize=(8,6))#,dpi=150)
	ax.plot(data2[:,per_col]./365.25 ,xprob(data2[:,end]),label="16 yrs")
	ax.plot(data4[:,per_col]./365.25 ,xprob(data4[:,end]),label="18 yrs")
	ax.plot(data6[:,per_col]./365.25 ,xprob(data6[:,end]),label="20 yrs")
	ax.plot(data8[:,per_col]./365.25 ,xprob(data8[:,end]),label="22 yrs")
	ax.plot(data10[:,per_col]./365.25 ,xprob(data10[:,end]),label="24 yrs")
	axvline(true_per,linestyle="--",color="black")
	text(true_per + true_per/100,1.01,pname,fontsize="large")
	#inset 
	# text(label_xloc,1.05,sim_obs_label)
	legend(loc="upper right",fontsize="large",title=case_label,title_fontsize="large",bbox_to_anchor=(0.,1.02,1.,.102),ncol=5,mode="expand",borderaxespad=0.0)
	if grid_type_nplanet=="p4"
	ax2=fig.add_axes([0.5,0.2,0.3,0.5])
	ax2.plot(data1[:,per_col]./365.25 ,xprob(data1[:,end]),label="15 yrs")
	ax2.plot(data2[:,per_col]./365.25 ,xprob(data2[:,end]),label="16 yrs")
	ax2.plot(data3[:,per_col]./365.25 ,xprob(data3[:,end]),label="17 yrs")
	ax2.plot(data4[:,per_col]./365.25 ,xprob(data4[:,end]),label="18 yrs")
	ax2.plot(data5[:,per_col]./365.25 ,xprob(data5[:,end]),label="19 yrs")
	ax2.plot(data6[:,per_col]./365.25 ,xprob(data6[:,end]),label="20 yrs")
	ax2.plot(data7[:,per_col]./365.25 ,xprob(data7[:,end]),label="21 yrs")
	ax2.plot(data8[:,per_col]./365.25 ,xprob(data8[:,end]),label="22 yrs")
	ax2.plot(data9[:,per_col]./365.25 ,xprob(data9[:,end]),label="23 yrs")
	ax2.plot(data10[:,per_col]./365.25 ,xprob(data10[:,end]),label="24 yrs")
	ax2.set_xlim(1.81,1.9)
	ax2.axvline(true_per,linestyle="--",color="black")
	ax2.tick_params(which="both",direction="in")
	end
	ax.set_xlabel("Planet Period Search Grid [years]",fontsize="x-large")
	# xlim(minimum(xgrid_in_yrs)-minimum(xgrid_in_yrs)/10,maximum(xgrid_in_yrs)+.05)
	ax.set_ylabel("Probability",fontsize="x-large")
	ax.set_ylim(0,1.2)
	ax.minorticks_on()
	tight_layout()
end
function moon_sigs(nyear,grid_type_nplanet,per_col,true_per,pname,case_num)
	file1=string("grid/",grid_type_nplanet,"_grid10s",nyear,"yrs.csv")
	file2=string("grid/",grid_type_nplanet,"_grid30s",nyear,"yrs.csv")
	file3=string("grid/",grid_type_nplanet,"_grid60s",nyear,"yrs.csv")
	file4=string("grid/",grid_type_nplanet,"_grid80s",nyear,"yrs.csv")
	file5=string("grid/",grid_type_nplanet,"_grid90s",nyear,"yrs.csv")
	file6=string("grid/",grid_type_nplanet,"_grid100s",nyear,"yrs.csv")
	file7=string("grid/",grid_type_nplanet,"_grid110s",nyear,"yrs.csv")
	file8=string("grid/",grid_type_nplanet,"_grid120s",nyear,"yrs.csv")
	case_label=string("Case 2"," [",nyear," yr span]")
	data1,header1=readdlm(file1,',',header=true)
	data2,header2=readdlm(file2,',',header=true)
	data3,header3=readdlm(file3,',',header=true)
	data4,header4=readdlm(file4,',',header=true)
	data5,header5=readdlm(file5,',',header=true)
	data6,header6=readdlm(file6,',',header=true)
	data7,header7=readdlm(file7,',',header=true)
	data8,header8=readdlm(file8,',',header=true)
fig,ax=subplots(figsize=(5,5))#,dpi=150)
ax.plot(data1[:,per_col] ,xprob(data1[:,end]),label=L"$\sigma_{obs}=$ 10 sec")
ax.plot(data2[:,per_col] ,xprob(data2[:,end]),label="30 sec")
ax.plot(data3[:,per_col] ,xprob(data3[:,end]),label="60 sec",linestyle="--")
ax.plot(data4[:,per_col] ,xprob(data4[:,end]),label="80 sec",linestyle="-.")
ax.plot(data5[:,per_col] ,xprob(data5[:,end]),label="90 sec",linestyle="-.")
ax.plot(data6[:,per_col] ,xprob(data6[:,end]),label="100 sec",linestyle="-.")
ax.plot(data7[:,per_col] ,xprob(data7[:,end]),label="110 sec",linestyle="-.")
ax.plot(data8[:,per_col] ,xprob(data8[:,end]),label="120 sec",linestyle="-.")
axvline(true_per,linestyle="--",color="black")
text(true_per + true_per/100,1.01,pname,fontsize="large")
# text(label_xloc,1.05,sim_obs_label)
legend(loc="upper right",fontsize="medium",title=case_label,title_fontsize="medium",bbox_to_anchor=(0.,1.02,1.,.102),ncol=4,mode="expand",borderaxespad=0.0)
xlabel("Satellite-induced Phase Offset Search Grid [rads]",fontsize="x-large")
# xlim(minimum(xgrid_in_yrs)-minimum(xgrid_in_yrs)/10,maximum(xgrid_in_yrs)+.05)
ylabel("Probability",fontsize="x-large")
ylim(0,1.2)
minorticks_on()
tick_params(which="both",direction="in")
tight_layout()
end
function moon_yrs(sigma,grid_type_nplanet,per_col,true_per,pname,case_num)
	file1=string("grid/",grid_type_nplanet,"_grid",sigma,"s15yrs.csv")
	file2=string("grid/",grid_type_nplanet,"_grid",sigma,"s16yrs.csv")
	file3=string("grid/",grid_type_nplanet,"_grid",sigma,"s17yrs.csv")
	file4=string("grid/",grid_type_nplanet,"_grid",sigma,"s18yrs.csv")
	file5=string("grid/",grid_type_nplanet,"_grid",sigma,"s19yrs.csv")
	file6=string("grid/",grid_type_nplanet,"_grid",sigma,"s20yrs.csv")
	file7=string("grid/",grid_type_nplanet,"_grid",sigma,"s21yrs.csv")
	file8=string("grid/",grid_type_nplanet,"_grid",sigma,"s22yrs.csv")
	file9=string("grid/",grid_type_nplanet,"_grid",sigma,"s23yrs.csv")
	file10=string("grid/",grid_type_nplanet,"_grid",sigma,"s24yrs.csv")
	# file10=string("grid/",grid_type_nplanet,"_grid",sigma,"s25yrs.csv")
	# file10=string("grid/",grid_type_nplanet,"_grid",sigma,"s26yrs.csv")
	case_label=string("Case 2"," [σ=",sigma," sec]")
	data1,header1=readdlm(file1,',',header=true)
	data2,header2=readdlm(file2,',',header=true)
	data3,header3=readdlm(file3,',',header=true)
	data4,header4=readdlm(file4,',',header=true)
	data5,header5=readdlm(file5,',',header=true)
	data6,header6=readdlm(file6,',',header=true)
	data7,header7=readdlm(file7,',',header=true)
	data8,header8=readdlm(file8,',',header=true)
	data9,header9=readdlm(file9,',',header=true)
	data10,header10=readdlm(file10,',',header=true)
fig,ax=subplots(figsize=(8,6))#,dpi=150)
ax.plot(data1[:,per_col] ,xprob(data1[:,end]),label="15 yrs")
# ax.plot(data2[:,per_col] ,xprob(data2[:,end]),label="16 yrs")
ax.plot(data3[:,per_col] ,xprob(data3[:,end]),label="17 yrs")
# ax.plot(data4[:,per_col] ,xprob(data4[:,end]),label="18 yrs")
ax.plot(data5[:,per_col] ,xprob(data5[:,end]),label="19 yrs")
# ax.plot(data6[:,per_col] ,xprob(data6[:,end]),label="20 yrs")
ax.plot(data7[:,per_col] ,xprob(data7[:,end]),label="21 yrs")
# ax.plot(data8[:,per_col] ,xprob(data8[:,end]),label="22 yrs")
ax.plot(data9[:,per_col] ,xprob(data9[:,end]),label="23 yrs")
# ax.plot(data10[:,per_col] ,xprob(data10[:,end]),label="24 yrs")
axvline(true_per,linestyle="--",color="black")
text(true_per + true_per/100,1.01,pname,fontsize="large")
#inset 
# text(label_xloc,1.05,sim_obs_label)
# legend(loc="upper right",fontsize="medium",title=case_label,title_fontsize="medium",bbox_to_anchor=(0.,1.02,1.,.102),ncol=5,mode="expand",borderaxespad=0.0)
legend(loc="upper right",fontsize="medium",title=case_label,title_fontsize="medium",bbox_to_anchor=(0.,0.9,1.,.102),ncol=5,mode="expand",borderaxespad=0.0)
xlabel("Satellite-induced Phase Offset Search Grid [rads]",fontsize="x-large")
# xlim(minimum(xgrid_in_yrs)-minimum(xgrid_in_yrs)/10,maximum(xgrid_in_yrs)+.05)
ylabel("Probability",fontsize="x-large")
ylim(0,1.2)
xlim(2.26,2.36)
minorticks_on()
tick_params(which="both",direction="in")
tight_layout()
end