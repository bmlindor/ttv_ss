using PyPlot,PyCall,Statistics,JLD2,DelimitedFiles,DataFrames
rc("font",family="sans-serif")
rc("lines",linewidth=2)
include("histogram.jl")

xprob(lprob)=exp.(lprob .- maximum(lprob))
 	# true_per3=1.8808476
 	# true_per4=11.862615
	# prob(p3,p3prob,p3_cur,nbins,param,"firebrick","Jupiter")
# Create likelihood/probability plots of search grids
function per_grid(sigma,nyear,grid_type_nplanet,per_col,true_per,pname,case_num,label_xloc)
	if case_num==1
	file=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 1"
	elseif case_num==2
	file=string("grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 2"
	end
	if grid_type_nplanet=="p3" || grid_type_nplanet=="widep3"
    model=L"$\mathcal{H}_{PPP}$"
    # per_col,true_per,pname,color=12,11.862615,"Jupiter","firebrick"; label_xloc=5
  elseif grid_type_nplanet=="p4" || grid_type_nplanet=="widep4"
    model=L"$\mathcal{H}_{PPPP}$"
    # per_col,true_per,pname,color=12,1.8808476,"Mars","orange"; label_xloc=2.75
  elseif grid_type_nplanet=="p3moon"
    model=L"$\mathcal{H}_{PPsP}$"
    # per_col,true_per,pname,color=12,11.862615,"Jupiter","firebrick"; label_xloc=5
  elseif grid_type_nplanet=="p3moonp4"
    model=L"$\mathcal{H}_{PPsPP}$"
    # per_col,true_per,pname,color=12,1.8808476,"Mars","orange"; label_xloc=2.75
  end
	fit,header=readdlm(file,',',header=true)
	sim_obs_label= string(case_label," [",nyear," yr span]")
	# save_as =string("IMAGES/wide_grids/case",case_num,"_",grid_type_nplanet,pname,sigma,"secs",nyear,"yrs.png")
	title(string(model," [",nyear," yr span]"))
	axvline(true_per,linestyle="--",color="black")
	text(1.1*true_per,1.01,pname,fontsize="large")
	#text(label_xloc,1.05,sim_obs_label,fontsize="medium")
 	plot(fit[:,per_col]./365.25,xprob(fit[:,end]),label=string(L"$\sigma_{obs}=$",sigma," sec"))#label=string("[",nyear," yr span]",L"$\sigma_{obs}=$",sigma," sec"))
 	xlabel("Planet Period Search Grid",fontsize="x-large")
 	ylabel("Probability",fontsize="x-large")
 	legend()
 	# ylim(0,1.2)
end
function moon_grid(sigma,nyear,grid_type_nplanet,per_col,true_per,color,pname,case_num,label_xloc)
	if case_num==1
	file=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 1"
	elseif case_num==2
	file=string("grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 2"
	end
	fit,header=readdlm(file,',',header=true)
	sim_obs_label= string(case_label," [",nyear," yr span]",'\n',L"$\sigma_{obs}=$",sigma," sec")
	# save_as =string("IMAGES/wide_grids/case",case_num,"_",grid_type_nplanet,pname,sigma,"secs",nyear,"yrs.png")
	axvline(true_per,linestyle="--",color="black")
	text(true_per + true_per/100,1.01,pname,fontsize="large")
	text(label_xloc,1.05,sim_obs_label,fontsize="medium")
 	plot(fit[:,per_col],xprob(fit[:,end]),color)
 	xlabel("Planet Period Search Grid",fontsize="x-large")
 	ylabel("Probability",fontsize="x-large")
 	ylim(0.98,1.2)
end

function wide_grid(sigma,nyear,grid_type_nplanet,case_num,nbins,include_moon::Bool=false)
	if case_num==1
	fitfile=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	widefit=string("grid/fromEMB/wide",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 1"
	title="Search from Venus + EMB TTVs"
	elseif case_num==2
  fitfile=string("grid/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	widefit=string("grid/wide",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 2"
	title="Search from Venus + Earth TTVs"
	end
  if grid_type_nplanet=="p3"
    model=L"$\mathcal{H}_{PPP}$"
    per_col,true_per,pname,color=12,11.862615,"Jupiter","firebrick"; label_xloc=5
    text_loc=11.9
  elseif grid_type_nplanet=="p4"
    model=L"$\mathcal{H}_{PPPP}$"
    per_col,true_per,pname,color=12,1.8808476,"Mars","orange"; label_xloc=2.75
    text_loc=1.85
  elseif grid_type_nplanet=="p3moon"
    model=L"$\mathcal{H}_{PPsP}$"
    per_col,true_per,pname,color=12,11.862615,"Jupiter","firebrick"; label_xloc=5
    text_loc=11.9
  elseif grid_type_nplanet=="p3moonp4"
    model=L"$\mathcal{H}_{PPsPP}$"
    per_col,true_per,pname,color=12,1.8808476,"Mars","orange"; label_xloc=2.75
    text_loc=0.9*true_per
  end
	fit,header=readdlm(fitfile,',',header=true)
	wide,header=readdlm(widefit,',',header=true)
	sim_obs_label= string(" [",nyear," yr span]",'\n',L"$\sigma_{obs}=$",sigma," sec")
	fig=figure(figsize=(5,5),dpi=150)
	ax1=subplot(211)
	ax1.set_title(string(title),fontsize="x-large")
	ax1.text(label_xloc,0.86,sim_obs_label,fontsize="small")
	ax1.axvline(true_per,linestyle="--",color="black")
 	ax1.plot(wide[:,per_col]./365.25,xprob(wide[:,end]),color="black")
	ax1.axvspan(minimum(fit[:,per_col]./365.25),maximum(fit[:,per_col]./365.25),color="lightblue",alpha=0.45)
	ax1.minorticks_on()
	ax1.tick_params(which="both",direction="in")
 	# ax1.plot(wide[:,per_col]./365.25,xprob(wide[:,end]),".")
 	ax2=subplot(212)
	ax2.axvline(true_per,linestyle="--",color="black")
	ax2.text(text_loc,0.96,pname,fontsize="large")
 	ax2.plot(fit[:,per_col]./365.25,xprob(fit[:,end]),color,label=string(model))
	# plot(fit[:,per_col]./365.25,xprob(fit[:,end]),".",label="pts")
	ax2.legend()
 	xlabel("Planet Period Search Grid",fontsize="x-large")
 	ylabel("Probability",fontsize="x-large")
 	ax2.minorticks_on()
	ax2.tick_params(which="both",direction="in")
 	# xlim(1.8,1.95)
    # xlim(11,12.5)
  title=string("IMAGES/likelihood/case",case_num,"_",grid_type_nplanet,"_",sigma,"s",nyear,"yrs.png")
  savefig(title)
end

function plot_grid(sigma,nyear,grid_type_nplanet,per_col,true_per,pname,case_num,label_xloc)
	if case_num==1
	file=string("grid_wide/fromEMB/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 1"

	elseif case_num==2
	file=string("grid_wide/",grid_type_nplanet,"_grid",sigma,"s",nyear,"yrs.csv")
	case_label="Case 2"
		# case_label=string("Case 1"," [",nyear," yr span]")
	end
	fit,header=readdlm(file,',',header=true)
	sim_obs_label= string(case_label," [",nyear," yr span]",'\n',L"$\sigma_{obs}=$",sigma," sec")
	axvline(true_per,linestyle="--",color="black")
	text(true_per + true_per/100,1.01,pname)
	text(label_xloc,1.05,sim_obs_label)
	label=string(sigma,"s ",nyear," yrs")
 	plot(fit[:,per_col]./365.25,xprob(fit[:,end]),label=label)
 	legend()
 	ylim(0,1.2)
end

#### plot_profile(sigma,nyear,grid_type_nplanet,per_col,true_per,pname,case_num)
function plot_profile(sigma::Real,nyear::Real,label_xloc::Real,case_num::Int,nbins=Int,include_moon=false,p3m=false)
if case_num==1
	fitfile3=string("grid/fromEMB/","p4","_grid",sigma,"s",nyear,"yrs.csv")
	widefit3=string("grid/fromEMB/wide","p4","_grid",sigma,"s",nyear,"yrs.csv")
  	mcfile3=string("MCMC/fromEMB/","p4","_mcmc",sigma,"s",nyear,"yrs.jld2")
	fitfile4=string("grid/fromEMB/","p3","_grid",sigma,"s",nyear,"yrs.csv")
	widefit4=string("grid/fromEMB/wide","p3","_grid",sigma,"s",nyear,"yrs.csv")
  	mcfile4=string("MCMC/fromEMB/","p3","_mcmc",sigma,"s",nyear,"yrs.jld2")
	case_label="Case 1"
	title="Search from Venus + EMB TTVs"
	label1=L"$\mathcal{H}_{PPPP}$ "
elseif case_num==2
case_label="Case 2"
  	if include_moon
  		fitfile3=string("grid/","p3moonp4","_grid",sigma,"s",nyear,"yrs.csv")
			widefit3=string("grid/wide","p3moonp4","_grid",sigma,"s",nyear,"yrs.csv")
  		mcfile3=string("MCMC/","p3moonp4","_mcmc",sigma,"s",nyear,"yrs.jld2")
  		title="Search from Venus + Earth TTVs"
		 	label1=L"$\mathcal{H}_{PPsPP}$ "
  	else
  		fitfile3=string("grid/","p4","_grid",sigma,"s",nyear,"yrs.csv")
			widefit3=string("grid/wide","p4","_grid",sigma,"s",nyear,"yrs.csv")
  		mcfile3=string("MCMC/","p4","_mcmc",sigma,"s",nyear,"yrs.jld2")
  		title="Search from Venus + Earth TTVs"
  		label1=L"$\mathcal{H}_{PPPP}$ "
  	end
if p3m
    fitfile4=string("grid/","p3moon","_grid",sigma,"s",nyear,"yrs.csv")
    widefit4=string("grid/wide","p3","_grid",sigma,"s",nyear,"yrs.csv")
    mcfile4=string("MCMC/","p3moon","_mcmc",sigma,"s",nyear,"yrs.jld2")
    println("Hppmp fit on right.")
else
 	
    fitfile4=string("grid/","p3","_grid",sigma,"s",nyear,"yrs.csv")
    widefit4=string("grid/wide","p3","_grid",sigma,"s",nyear,"yrs.csv")
    mcfile4=string("MCMC/","p3","_mcmc",sigma,"s",nyear,"yrs.jld2")
end

end
	mc3=jldopen(String(mcfile3),"r")
 	mc4=jldopen(String(mcfile4),"r")
	fit3,header=readdlm(fitfile3,',',header=true)
	wide3,header=readdlm(widefit3,',',header=true)
	fit4,header=readdlm(fitfile4,',',header=true)
	wide4,header=readdlm(widefit4,',',header=true)

	par_mcmc3=vec(mc3["par_mcmc"][:,mc3["iburn"]:end,12])./365.25
	par_mcmc4=vec(mc4["par_mcmc"][:,mc4["iburn"]:end,12])./365.25
	sigsys3=(median(vec(mc3["par_mcmc"][:,mc3["iburn"]:end,end]))).* 3600*24
	sigsys_err3=(std(vec(mc3["par_mcmc"][:,mc3["iburn"]:end,end]))).* 3600*24
  sigtot3=sqrt(sigsys3^2 + sigma^2) 
	sigsys4=(median(vec(mc4["par_mcmc"][:,mc4["iburn"]:end,end]))).* 3600*24
	sigsys_err4=(std(vec(mc4["par_mcmc"][:,mc4["iburn"]:end,end]))).* 3600*24
	sigtot4=sqrt(sigsys4^2 + sigma^2) 

 	true_per3=1.8808476
 	true_per4=11.862615
	sim_obs_label= string(nyear," yr span",'\n',L"$\sigma_{obs}=$",sigma,"s")
	fig=figure(figsize=(5,5),dpi=150)
	# fig=figure(figsize=(5,5))
	# subplots_adjust(hspace=0.2,wspace=0.15)
	# ax1.title()
	ax1=subplot(211)
	ax1.set_title(string(title),fontsize="x-large")
	axvline(true_per3,linestyle="--",color="black")
 	axvline(true_per4,linestyle="--",color="black")
	text(label_xloc,0.8,sim_obs_label,fontsize="medium")
	# text(true_per + true_per/100,1.01,pname,fontsize="large")
	# text(minimum(wide[:,per_col]./365.25),1.05,sim_obs_label,fontsize="medium")
	ax1.axvspan(1.8,2,color="lightblue",alpha=0.45)
	ax1.axvspan(11.0,12.5,color="lightblue",alpha=0.45)
 	ax1.plot(wide3[:,12]./365.25,xprob(wide3[:,end]),color="orange")
	ax1.plot(wide4[:,12]./365.25,xprob(wide4[:,end]),color="firebrick")
 	ax1.set_xlabel("Planet Period [yrs]",fontsize="x-large")
 	ax1.set_ylabel("Probability",fontsize="x-large")
 	# ylim(0,1.15)
 	minorticks_on()
	tick_params(which="both",direction="in")

	## plot Mars zoom-in
	if include_moon
	ax2=fig.add_subplot(223,title=string(label1,L"$\sigma_{sys}=$",round(sigsys3,sigdigits=4),"s"))
	title=string("IMAGES/likelihood/wide_case",case_num,"_",sigma,"s",nyear,".png")
	else
	ax2=fig.add_subplot(223,title=string(L"$\mathcal{H}_{PPPP}$ ",L"$\sigma_{sys}=$",round(sigsys3,sigdigits=4)," s"))
	title=string("IMAGES/likelihood/wide_case",case_num,"_",sigma,"s",nyear,"_p4.png")
	end
	# text(1.83,1.1,"a)",fontweight="bold")
	axvline(true_per3,linestyle="--",color="black")
	text(1.89,0.96,"Mars",fontsize="medium")
	xbin,xhist,xbin_square,hist_square=histogram(par_mcmc3,nbins)
	# ax2.hist(par_mcmc3,bins=nbins,histtype="step",density=tr,color="orange")
 	ax2.plot(fit3[:,12]./365.25,xprob(fit3[:,end]),color="orange",label="Fit")
 	ax2.plot(xbin_square,hist_square./maximum(hist_square),color="orange",label="Posterior",alpha=0.75)
 	minorticks_on()
 	# legend()
	tick_params(which="both",direction="in")
	ax2.set_xlim(1.83,1.95)
		# ax2.set_xlim(quantile(par_mcmc3,0.1587))
	ax2.set_xlabel(L"Per$_4$ [yrs]",fontsize="large")
	ylabel("Probability",fontsize="x-large")

 	## plot Jupiter zoom-in
 	ax3=fig.add_subplot(224,sharey=ax2,title=string(L"$\mathcal{H}_{PPP}$ ",L"$\sigma_{sys}=$",round(sigsys4,sigdigits=4),"s"))
	# text(11,1.1,"b)",fontweight="bold")
 	axvline(true_per4,linestyle="--",color="black")
	text(12,0.96,"Jupiter",fontsize="medium")
	xbin,xhist,xbin_square,hist_square=histogram(par_mcmc4,50)
	# ax3.hist(par_mcmc4,bins=nbins,histtype="step",density=tr,color="firebrick")
 	ax3.plot(fit4[:,12]./365.25,xprob(fit4[:,end]),color="firebrick",label="Fit")
 	ax3.plot(xbin_square,hist_square./maximum(hist_square),color="firebrick",label="Poster",alpha=0.75)
 	minorticks_on()
	# legend(loc="upper right",fontsize="medium",title="Jupiter",title_fontsize="medium",bbox_to_anchor=(0.,1.02,1.,.102),ncol=5,mode="expand",borderaxespad=0.0)	
	ax3.set_xlim(11,12.5)
	ax3.tick_params(which="both",direction="in",left=true,labelleft=false)
	ax3.set_xlabel(L"Per$_3$ [yrs]",fontsize="large")
	# ax3.legend([m4])
	tight_layout()
 	savefig(title)
	# if include_moon
	# 	phi_col,true_phi,color=18,2.31586,"purple"; label_xloc=2.75
	# 	# lim=f["dpin"],f["dpout"]
	# 	fig=plt.figure()
	#   ax3=subplot(211)
	# 	ax3.axvline(true_phi,linestyle="--",color="black")
	# 	ax3.text(2.4,1,"Moon",fontsize="large")
	# 	ax3.plot(wide3[:,23]./365.25,xprob(wide3[:,end]),color="black")
	# 	# tight_layout()
	# 	# ax3.plot(grid_wide,lprob_wide,color="red") 
	# 	xlabel("Satellite Phase Offset Search Grid [radians]")
	# 	ylabel("Probability")
	# 	ax3.minorticks_on()
	# 	ax3.tick_params(which="both",direction="in")
	# 	# Inset zoom of finer δϕ grid 
	# 	ax4=fig.add_axes([0.68,0.4,0.2,0.4])
	#   ax4.plot(fit[:,23]./365.25,xprob(fit[:,end]),color="purple",label=string(model)) 
	#   ax4.plot(xbin_square,hist_square./maximum(hist_square),color="purple",alpha=0.5)
	#   xlim(2.24,2.37)
	#   ax4.minorticks_on()
	#   ax4.tick_params(which="both",direction="in",left="false",labelleft="false",
	#   		right="true",labelright="true")
	#   #title2=string("IMAGES/likelihood/",grid_type_nplanet,"Moon",sigma,"s",nyear,"yrs.png")
	#   #savefig(title2)
	# end
	# close()
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
	fit1,header1=readdlm(file1,',',header=true)
	fit2,header2=readdlm(file2,',',header=true)
	fit3,header3=readdlm(file3,',',header=true)
	fit4,header4=readdlm(file4,',',header=true)
		# save_as =string("IMAGES/wide_grids/case",case_num,"_",grid_type_nplanet,pname,sigma,"secs",nyear,"yrs.png")
	fig,ax=subplots(figsize=(5,5))#,dpi=150)
	ax.plot(fit1[:,per_col]./365.25 ,xprob(fit1[:,end]),label=L"$\sigma_{obs}=$ 10 sec")
	ax.plot(fit2[:,per_col]./365.25 ,xprob(fit2[:,end]),label="30 sec")
	ax.plot(fit3[:,per_col]./365.25 ,xprob(fit3[:,end]),label="60 sec",linestyle="--")
	ax.plot(fit4[:,per_col]./365.25 ,xprob(fit4[:,end]),label="80 sec",linestyle="-.")
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
### compare_yrs(30,"p5",22,29.447,"Saturn",1)

function compare_yrs(sigma,grid_type_nplanet,case_num)
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
	file11=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s25yrs.csv")
	file12=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s26yrs.csv")
	file13=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s27yrs.csv")
	file14=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s28yrs.csv")
	file15=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s29yrs.csv")
	file16=string("grid/fromEMB/",grid_type_nplanet,"_grid",sigma,"s30yrs.csv")
	case_label=string("σ=",sigma,"s")
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
	file11=string("grid/",grid_type_nplanet,"_grid",sigma,"s25yrs.csv")
	file12=string("grid/",grid_type_nplanet,"_grid",sigma,"s26yrs.csv")
	file13=string("grid/",grid_type_nplanet,"_grid",sigma,"s27yrs.csv")
	file14=string("grid/",grid_type_nplanet,"_grid",sigma,"s28yrs.csv")
	file15=string("grid/",grid_type_nplanet,"_grid",sigma,"s29yrs.csv")
	file16=string("grid/",grid_type_nplanet,"_grid",sigma,"s30yrs.csv")
	case_label=string("σ=",sigma,"s")
	end
	# fit1,header1=readdlm(file1,',',header=true)
			if grid_type_nplanet=="p3" || grid_type_nplanet=="widep3"
	fit2,header2=readdlm(file2,',',header=true)
	fit3,header3=readdlm(file3,',',header=true)
	fit4,header4=readdlm(file4,',',header=true)
	fit5,header5=readdlm(file5,',',header=true)
	fit6,header6=readdlm(file6,',',header=true)
	fit7,header7=readdlm(file7,',',header=true)
		end
	fit8,header8=readdlm(file8,',',header=true)
	fit9,header9=readdlm(file9,',',header=true)
	fit10,header10=readdlm(file10,',',header=true)
	fit11,header11=readdlm(file11,',',header=true)
	fit12,header12=readdlm(file12,',',header=true)
	fit13,header13=readdlm(file13,',',header=true)
	fit14,header14=readdlm(file14,',',header=true)
	fit15,header15=readdlm(file15,',',header=true)
	fit16,header16=readdlm(file16,',',header=true)
	if grid_type_nplanet=="p3" || grid_type_nplanet=="widep3"
    model=L"$\mathcal{H}_{PPP}$"
    per_col,true_per,pname,color=12,11.862615,"Jupiter","firebrick";
    pnum="3rd"
  elseif grid_type_nplanet=="p4" || grid_type_nplanet=="widep4"
    model=L"$\mathcal{H}_{PPPP}$"
    per_col,true_per,pname,color=12,1.8808476,"Mars","orange";
    pnum="4th"
  elseif grid_type_nplanet=="p3moon"
    model=L"$\mathcal{H}_{PPsP}$"
    per_col,true_per,pname,color=12,11.862615,"Jupiter","firebrick";
    pnum="3rd"
  elseif grid_type_nplanet=="p3moonp4"
    model=L"$\mathcal{H}_{PPsPP}$"
    per_col,true_per,pname,color=12,1.8808476,"Mars","orange";
    pnum="4th"
  end

	fig,ax=subplots(figsize=(8,6))#,dpi=150)
	title(string("Search for ",pnum," Planet (i.e. ",pname,") in ",model),fontsize="xx-large")
	# ax.plot(fit1[:,per_col]./365.25 ,xprob(fit1[:,end]),label="15 yrs")
	if grid_type_nplanet=="p3" || grid_type_nplanet=="widep3"
	ax.plot(fit2[:,per_col]./365.25 ,xprob(fit2[:,end]),label="16 yrs",linewidth=1)
	ax.plot(fit3[:,per_col]./365.25 ,xprob(fit3[:,end]),label="17 yrs",linewidth=1)
	ax.plot(fit4[:,per_col]./365.25 ,xprob(fit4[:,end]),label="18 yrs",linewidth=1)
	ax.plot(fit5[:,per_col]./365.25 ,xprob(fit5[:,end]),label="19 yrs",linewidth=1)
	ax.plot(fit6[:,per_col]./365.25 ,xprob(fit6[:,end]),label="20 yrs",linewidth=1)
	ax.plot(fit7[:,per_col]./365.25 ,xprob(fit7[:,end]),label="21 yrs",linewidth=1)
	end
	ax.plot(fit8[:,per_col]./365.25 ,xprob(fit8[:,end]),label="22 yrs",linewidth=1)#,linestyle="--")	
	ax.plot(fit9[:,per_col]./365.25 ,xprob(fit9[:,end]),label="23 yrs",linewidth=1)#,linestyle="--")
	ax.plot(fit10[:,per_col]./365.25 ,xprob(fit10[:,end]),label="24 yrs",linewidth=1)#,linestyle="--")
	ax.plot(fit11[:,per_col]./365.25 ,xprob(fit11[:,end]),label="25 yrs",linewidth=1)#,linestyle="--")
	ax.plot(fit12[:,per_col]./365.25 ,xprob(fit12[:,end]),label="26 yrs",linewidth=1,linestyle="--")
	ax.plot(fit13[:,per_col]./365.25 ,xprob(fit13[:,end]),label="27 yrs",linewidth=1,linestyle="--")  
	ax.plot(fit14[:,per_col]./365.25 ,xprob(fit14[:,end]),label="28 yrs",linewidth=1,linestyle="--")  
	ax.plot(fit15[:,per_col]./365.25 ,xprob(fit15[:,end]),label="29 yrs",linewidth=1,linestyle="--")  
	ax.plot(fit16[:,per_col]./365.25 ,xprob(fit16[:,end]),label="30 yrs",linewidth=1,linestyle="--")      
	axvline(true_per,linestyle="--",color="black")
	text(1.01*true_per,0.96,pname,fontsize="large")
	tick_params(which="both",direction="in")
	# text(true_per + true_per/100,1.01,pname,fontsize="large")
	#inset 
	# text(label_xloc,1.05,sim_obs_label)
	ax.legend(loc="upper right",fontsize="large",title=case_label,title_fontsize="large",bbox_to_anchor=(0.,0.9,1.,.102),ncol=5,mode="expand",borderaxespad=0.0)
	# legend(fontsize="medium",title_fontsize="medium",title=case_label)
	# if grid_type_nplanet=="p4" #|| grid_type_nplanet=="p3moonp4"
	# ax2=fig.add_axes([0.5,0.2,0.3,0.5])
	# ax2.plot(fit1[:,per_col]./365.25 ,xprob(fit1[:,end]),label="15 yrs")
	# # ax2.plot(fit2[:,per_col]./365.25 ,xprob(fit2[:,end]),label="16 yrs")
	# ax2.plot(fit3[:,per_col]./365.25 ,xprob(fit3[:,end]),label="17 yrs")
	# # ax2.plot(fit4[:,per_col]./365.25 ,xprob(fit4[:,end]),label="18 yrs")
	# ax2.plot(fit5[:,per_col]./365.25 ,xprob(fit5[:,end]),label="19 yrs")
	# # ax2.plot(fit6[:,per_col]./365.25 ,xprob(fit6[:,end]),label="20 yrs")
	# ax2.plot(fit7[:,per_col]./365.25 ,xprob(fit7[:,end]),label="21 yrs")
	# # ax2.plot(fit8[:,per_col]./365.25 ,xprob(fit8[:,end]),label="22 yrs")
	# ax2.plot(fit9[:,per_col]./365.25 ,xprob(fit9[:,end]),label="23 yrs")
	# # ax2.plot(fit10[:,per_col]./365.25 ,xprob(fit10[:,end]),label="24 yrs")
	# ax2.set_xlim(1.81,1.9)
	# ax2.axvline(true_per,linestyle="--",color="black")
	# ax2.tick_params(which="both",direction="in")
	# end
	ax.set_xlabel("Planet Period Search Grid [years]",fontsize="xx-large")
	# # xlim(minimum(xgrid_in_yrs)-minimum(xgrid_in_yrs)/10,maximum(xgrid_in_yrs)+.05)
	ax.set_ylabel("Probability",fontsize="xx-large")
	ax.set_ylim(0,1.3)
	# #ax.set_xlim(1.8,1.9)
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
	fit2,header2=readdlm(file2,',',header=true)
	fit4,header4=readdlm(file4,',',header=true)
	fit6,header6=readdlm(file6,',',header=true)
	fit8,header8=readdlm(file8,',',header=true)
	fit10,header10=readdlm(file10,',',header=true)
	fig,ax=subplots(figsize=(8,6))#,dpi=150)
	ax.plot(fit2[:,per_col]./365.25 ,xprob(fit2[:,end]),label="16 yrs")
	ax.plot(fit4[:,per_col]./365.25 ,xprob(fit4[:,end]),label="18 yrs")
	ax.plot(fit6[:,per_col]./365.25 ,xprob(fit6[:,end]),label="20 yrs")
	ax.plot(fit8[:,per_col]./365.25 ,xprob(fit8[:,end]),label="22 yrs")
	ax.plot(fit10[:,per_col]./365.25 ,xprob(fit10[:,end]),label="24 yrs")
	axvline(true_per,linestyle="--",color="black")
	text(true_per + true_per/100,1.01,pname,fontsize="large")
	#inset 
	# text(label_xloc,1.05,sim_obs_label)
	legend(loc="upper right",fontsize="large",title=case_label,title_fontsize="large",bbox_to_anchor=(0.,1.02,1.,.102),ncol=5,mode="expand",borderaxespad=0.0)
	if grid_type_nplanet=="p4"
	ax2=fig.add_axes([0.5,0.2,0.3,0.5])
	ax2.plot(fit1[:,per_col]./365.25 ,xprob(fit1[:,end]),label="15 yrs")
	ax2.plot(fit2[:,per_col]./365.25 ,xprob(fit2[:,end]),label="16 yrs")
	ax2.plot(fit3[:,per_col]./365.25 ,xprob(fit3[:,end]),label="17 yrs")
	ax2.plot(fit4[:,per_col]./365.25 ,xprob(fit4[:,end]),label="18 yrs")
	ax2.plot(fit5[:,per_col]./365.25 ,xprob(fit5[:,end]),label="19 yrs")
	ax2.plot(fit6[:,per_col]./365.25 ,xprob(fit6[:,end]),label="20 yrs")
	ax2.plot(fit7[:,per_col]./365.25 ,xprob(fit7[:,end]),label="21 yrs")
	ax2.plot(fit8[:,per_col]./365.25 ,xprob(fit8[:,end]),label="22 yrs")
	ax2.plot(fit9[:,per_col]./365.25 ,xprob(fit9[:,end]),label="23 yrs")
	ax2.plot(fit10[:,per_col]./365.25 ,xprob(fit10[:,end]),label="24 yrs")
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
	fit1,header1=readdlm(file1,',',header=true)
	fit2,header2=readdlm(file2,',',header=true)
	fit3,header3=readdlm(file3,',',header=true)
	fit4,header4=readdlm(file4,',',header=true)
	fit5,header5=readdlm(file5,',',header=true)
	fit6,header6=readdlm(file6,',',header=true)
	fit7,header7=readdlm(file7,',',header=true)
	fit8,header8=readdlm(file8,',',header=true)
	fig,ax=subplots(figsize=(5,5))#,dpi=150)
	ax.plot(fit1[:,per_col] ,xprob(fit1[:,end]),label=L"$\sigma_{obs}=$ 10 sec")
	ax.plot(fit2[:,per_col] ,xprob(fit2[:,end]),label="30 sec")
	ax.plot(fit3[:,per_col] ,xprob(fit3[:,end]),label="60 sec",linestyle="--")
	ax.plot(fit4[:,per_col] ,xprob(fit4[:,end]),label="80 sec",linestyle="-.")
	ax.plot(fit5[:,per_col] ,xprob(fit5[:,end]),label="90 sec",linestyle="-.")
	ax.plot(fit6[:,per_col] ,xprob(fit6[:,end]),label="100 sec",linestyle="-.")
	ax.plot(fit7[:,per_col] ,xprob(fit7[:,end]),label="110 sec",linestyle="-.")
	ax.plot(fit8[:,per_col] ,xprob(fit8[:,end]),label="120 sec",linestyle="-.")
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
	fit1,header1=readdlm(file1,',',header=true)
	fit2,header2=readdlm(file2,',',header=true)
	fit3,header3=readdlm(file3,',',header=true)
	fit4,header4=readdlm(file4,',',header=true)
	fit5,header5=readdlm(file5,',',header=true)
	fit6,header6=readdlm(file6,',',header=true)
	fit7,header7=readdlm(file7,',',header=true)
	fit8,header8=readdlm(file8,',',header=true)
	fit9,header9=readdlm(file9,',',header=true)
	fit10,header10=readdlm(file10,',',header=true)
	fig,ax=subplots(figsize=(8,6))#,dpi=150)
	ax.plot(fit1[:,per_col] ,xprob(fit1[:,end]),label="15 yrs")
	# ax.plot(fit2[:,per_col] ,xprob(fit2[:,end]),label="16 yrs")
	ax.plot(fit3[:,per_col] ,xprob(fit3[:,end]),label="17 yrs")
	# ax.plot(fit4[:,per_col] ,xprob(fit4[:,end]),label="18 yrs")
	ax.plot(fit5[:,per_col] ,xprob(fit5[:,end]),label="19 yrs")
	# ax.plot(fit6[:,per_col] ,xprob(fit6[:,end]),label="20 yrs")
	ax.plot(fit7[:,per_col] ,xprob(fit7[:,end]),label="21 yrs")
	# ax.plot(fit8[:,per_col] ,xprob(fit8[:,end]),label="22 yrs")
	ax.plot(fit9[:,per_col] ,xprob(fit9[:,end]),label="23 yrs")
	# ax.plot(fit10[:,per_col] ,xprob(fit10[:,end]),label="24 yrs")
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
