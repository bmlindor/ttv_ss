
using PyPlot
function histogram(param,nbin)
  p1 = minimum(param)-1e-15; p2 = maximum(param)+1e-15
  pbin_square = [p1]
  hist_square = [0.0]
  pbin = zeros(nbin)
  hist = zeros(nbin)
  psort = sort(param)
  i1 = 1; np = size(param)[1]
  #println(p1," ",psort[i1]," ",p2," ",psort[np])
  for i=1:nbin
    while psort[i1] <= p1+(p2-p1)/nbin*i
      hist[i] += 1.0
      i1 += 1
      if i1 == np
        break
      end
    end
    push!(hist_square,hist[i]); push!(hist_square,hist[i])
    pbin[i] = p1+(p2-p1)/nbin*(i-0.5)
    push!(pbin_square,p1+(p2-p1)/nbin*(i-1))
    push!(pbin_square,p1+(p2-p1)/nbin*i)
    if i1 == np
      break
    end
  end
  push!(hist_square,0.0)
  push!(pbin_square,p2)
  return pbin,hist,pbin_square,hist_square
end

function filled_hist(ax,param,nbin,label,color::String)
	# ax.fill_between(xbin_square, xhist_square,alpha=0.25,label=label,color=color)
  ax.hist(param,nbin,histtype="step",fill=true,facecolor=color,label=label,edgecolor=color,density=true,linewidth=2,zorder=1)
	return
end
function lined_hist_stack(ax,stacked_data,labels=nothing,linestyles=nothing,colors=nothing;nbins=50,line::Bool=true)
  if isnothing(labels)
  labels=["label $i" for i=1:length(data)]
  end
  if isnothing(colors)
    julia_colors=["#389826","#CB3C33","#9558B2","#4063D8"]
    colors=["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b"]
  end
  if isnothing(linestyles)
    linestyles=["--","-","-."]
  end
  for (j,(data,label,ls,color)) in enumerate(zip(stacked_data,labels,linestyles,colors))
  if label==""
  continue
  end
    ax.hist(data,nbins,histtype="step",alpha=0.4,density=true,color=color,linewidth=1.5,label=label,linestyle=ls,fill=true,facecolor=color)
  end
end
function stacked_hist(ax,stacked_data,labels=nothing,colors=nothing;nbins=50,fill::Bool=true)
  if isnothing(labels)
    labels=["label $i" for i=1:length(data)]
  end

  if isnothing(colors)
    julia_colors=["#389826","#CB3C33","#9558B2","#4063D8"]
    colors=["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b"]
  end

  for (j,(data,label,color)) in enumerate(zip(stacked_data,labels,colors))
  if label==""
  continue
  end
    # ax.hist(data,10,histtype="step",alpha=1.0,density=true,zorder=1)
    filled_hist(ax,data,nbins,label,color)
  end
  # ax.legend(fontsize="medium")
end

function comp_hist(sigma,nyear,grid_type_nplanet,nbins,case=1,include_moon=false;grid_type_nplanet2="p2",grid_type_nplanet3="p3",grid_type_nplanet4="p4")
  if  grid_type_nplanet=="p2" 
      model=L"$\mathcal{H}_{PP}$"
    elseif grid_type_nplanet=="p3" || grid_type_nplanet=="widep3"
      model=L"$\mathcal{H}_{PPP}$"
    elseif grid_type_nplanet=="p4" || grid_type_nplanet=="widep4"
      model=L"$\mathcal{H}_{PPPP}$"
    elseif grid_type_nplanet=="p3moon" || grid_type_nplanet=="widep3moon"
      model=L"$\mathcal{H}_{PPsP}$"
    elseif grid_type_nplanet=="p3moonp4" || grid_type_nplanet=="widep3moonp4"
      model=L"$\mathcal{H}_{PPsPP}$"
  end  
   parname=[
    L"$m_b / M_{⋆}$",L"$P_b$ [days]",L"$t_{0,b}$ [JD $- 2.43e6$]",L"$e_b cos(ω_b)$",L"$e_b sin(ω_b)$",
    L"$m_c / M_{⋆}$",L"$P_c$ [days]",L"$t_{0,c}$ [JD $- 2.43e6$]",L"$e_c cos(ω_c)$",L"$e_c sin(ω_c)$",
    L"$m_e / M_{⋆}$",L"$P_e$ [days]",L"$t_{0,e}$ [JD $- 2.43e6$]",L"$e_e cos(ω_e)$",L"$e_e sin(ω_e)$",
    L"$m_d / M_{⋆}$",L"$P_d$ [days]",L"$t_{0,d}$ [JD $- 2.43e6$]",L"$e_d cos(ω_d)$",L"$e_d sin(ω_d)$",
    L"$μ_f$",L"$P_f$ [days]",L"$t_{0,f}$",L"$e_f cos(ω_f)$",L"$e_f sin(ω_f)$"]
   model4=L"$\mathcal{H}_{PPPP}$"
    model2=L"$\mathcal{H}_{PP}$"
  model3=L"$\mathcal{H}_{PPP}$"
   mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  mcfile2=string("MCMC/fromEMB/",grid_type_nplanet2,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  mcfile3=string("MCMC/fromEMB/",grid_type_nplanet3,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  function make_plot(ax,param_col,label,linestyle=nothing,color=nothing;nbins=50)
    values=[];labels=[]
    if isfile(mcfile)
    jldmc=jldopen(mcfile,"r");    jldmc2=jldopen(mcfile2,"r");    jldmc3=jldopen(mcfile3,"r")
    # nwalkers,nsteps=jldmc["nwalkers"],jldmc["nsteps"]
    # iburn,samples=jldmc["iburn"], jldmc["indepsamples"] param=jldmc["param"]
    # pname=jldmc["pname"]
    if param_col <= 10
    values=[ vec(jldmc2["par_mcmc"][:,jldmc2["iburn"]:end,param_col]),
    vec(jldmc3["par_mcmc"][:,jldmc3["iburn"]:end,param_col]), 
    vec(jldmc["par_mcmc"][:,jldmc["iburn"]:end,param_col])]
    labels=["$model2","$model3","$model4"]
    elseif 11 <= param_col <= 15 # plot jup
    values=[
    vec(jldmc3["par_mcmc"][:,jldmc3["iburn"]:end,param_col]),   
    vec(jldmc["par_mcmc"][:,jldmc["iburn"]:end,param_col+5])]
    labels=["$model3","$model4"]
    elseif 16<= param_col <= 20 # plot mars
    values=[ 
    vec(jldmc["par_mcmc"][:,jldmc["iburn"]:end,param_col-5])]
    labels=["$model4"]
    end
    ax.ticklabel_format(useMathText=true)
    # ax.tick_params(bottom=false,top=true,labeltop=true,labelbottom=false)
    ax.set_title(label,loc="right")
    return lined_hist_stack(ax,values,labels,linestyle,color;nbins)
    end
  end
  fig,axs=subplots(3,5,figsize=(12,6.5))
  make_plot(axs[1,1],1,parname[1])
  fig.legend(title="Models",title_fontsize="x-large",loc="upper left",fontsize="large",bbox_to_anchor=(0.0,0.88,0.4,0.1))
  make_plot(axs[1,2],2,parname[2])
  make_plot(axs[1,3],3,parname[3]) # t0
  make_plot(axs[1,4],4,parname[4])
  make_plot(axs[1,5],5,parname[5])
  make_plot(axs[2,1],6,parname[6])
  make_plot(axs[2,2],7,parname[7])
  make_plot(axs[2,3],8,parname[8]) # t0
  make_plot(axs[2,4],9,parname[9])
  make_plot(axs[2,5],10,parname[10])
  make_plot(axs[3,1],11,parname[16],["-","-."],["#ff7f0e","#2ca02c"])
  make_plot(axs[3,2],12,parname[17],["-","-."],["#ff7f0e","#2ca02c"])
  make_plot(axs[3,3],13,parname[18],["-","-."],["#ff7f0e","#2ca02c"]) # t0
  make_plot(axs[3,4],14,parname[19],["-","-."],["#ff7f0e","#2ca02c"])
  make_plot(axs[3,5],15,parname[20],["--","-."],["#ff7f0e","#2ca02c"])

  #  # make_plot(axs[4,1],16,parname[11],["-."],["#2ca02c"])
   # make_plot(axs[4,2],17,parname[12],["-."],["#2ca02c"];nbins=1000)
   # make_plot(axs[4,3],19,parname[14],["-."],["#2ca02c"])
   # make_plot(axs[4,4],20,parname[15],["-."],["#2ca02c"])
   # axs[4,2].set_xlim(650,800);   axs[4,4].set_xlim(-0.25,0.05);
  axs[1,5].set_xlim(-0.05,0.05);  axs[2,5].set_xlim(-0.025,0.025);   axs[3,5].set_xlim(-0.1,0.025);
  axs[1,4].set_xlim(-0.05,0.05); axs[2,4].set_xlim(-0.05,0.05);   axs[3,4].set_xlim(-0.1,0.1);
  axs[2,1].set_xlim(2.35e-6,3.5e-6);  axs[1,1].set_xlim(1.6e-6,3.5e-6);
  # fig.suptitle
  fig.subplots_adjust(wspace=0.4,hspace=0.5,bottom=0.05,right=0.98,top=0.92)
   # tight_layout()
   title=string("IMAGES/discussion/case",case,"_",grid_type_nplanet,"_",sigma,"s",nyear,"yrs_common1D.png")
  savefig(title,dpi=200)
end


# cmap=plt.cm.get_cmap("plasma")
# ls_cycle=plt.cycler("linestyle",["-","--","-"])
# sty_cycle=col_cycle+ls_cycle
# fig,ax=subplots(2)
# N=1000;x = randn(N);y=randn(N);z=randn(N)
# data=[x,y,z]
# col_cycle=plt.cycler("color",cmap(range(0,1,length=length(data))))
# # ax[1].set_prop_cycle(ls_cycle)
# # plt.rc('axes', prop_cycle=default_cycler)
# stacked_hist(ax[1],data,["a","b","c"],["salmon","forestgreen","firebrick"])
# ax[2].hist(x, 10, histtype="step",fill=true,facecolor="lightgreen",label="t1",edgecolor="green",alpha=0.4)
# ax[2].hist(y, 10, histtype="step",fill=true,facecolor="lightgreen",label="t2",edgecolor="green",alpha=0.4)
# ax[2].legend()
