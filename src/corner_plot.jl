using JLD2,PyPlot,Statistics,Distributions,LinearAlgebra,PyCall
#include("CGS.jl")
include("MCMC.jl")
np = pyimport("numpy")
@pyimport matplotlib.colors as pltco
# rc("font",family="sans-serif")
# rc("lines",linewidth=2)
        
# Corner plotting method, similar to corner.py
function corner(xs,labels,bins=5;quantiles)
  K=length(xs)
  factor = 1.5  # size of one side of one panel
  lbdim = 0.5 * factor  # size of left/bottom margin
  trdim = 0.2 * factor  # size of top/right margin
  whspace = 0.05 * 2  # w/hspace size
  plotdim = factor * K + factor * (K - 1.0) * whspace
  dim = lbdim + plotdim + trdim
  # @show dim

  fig, axes = plt.subplots(K, K, figsize=(dim, dim))
  lb = lbdim / dim
  tr = (lbdim + plotdim) / dim
  fig.subplots_adjust(
      left=lb, bottom=lb, right=tr, top=tr, wspace=whspace, hspace=whspace)

  for (i ,x) in enumerate(xs)
      ax=axes[i,i]
      ax.hist(x,bins=bins,histtype="step",linewidth=2)
  #     ax.tick_params(left=false,labelleft=false)
      if length(quantiles)>0
          qvalues=quantile(x,quantiles)
          ax.axvspan(qvalues[1],qvalues[2],linestyle="--",color="limegreen",alpha=0.3)            
          for q in qvalues
              ax.axvline(q,linestyle="--",color="k")
          end
      end
      q_mid = median(x) 
      q_lo,q_hi=quantile(x,quantiles)
      q_m, q_p = q_mid - q_lo, q_hi - q_mid
  #     fmt = "{{0:{0}}}".format(title_fmt).format
  #     title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
  #     title = title.format(fmt(q_mid), fmt(q_m), fmt(q_p))
  #     title=string(L"$=$",q_mid,L"$_{}")
  #     ax.set_title(title)
       ax.set_xticklabels([])
       ax.set_yticklabels([])
      ax.set_yticks([])

      for (j, y) in enumerate(xs)
          ax=axes[i,j]
           if j > i
                  ax.set_frame_on(false)
                  ax.set_xticks([])
                  ax.set_yticks([])
                  continue
              elseif j == i
                  continue
              end
          h,xedges,yedges=ax.hist2d(y,x,bins=bins,alpha=0.9)
          h_sort=sort(vec(h))
        sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
        sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
      level1=h_sort[sig2[1]];level2=h_sort[sig1[1]]
      if level1==level2
      error("Too many bins or not enough data points")
      end
      ax.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[level1,level2])
      if length(labels)>0
          if i < K 
          ax.set_xticklabels([])
          else
          ax.set_xlabel(labels[j],fontsize="medium")
          end
          if j>1
          ax.set_yticklabels([])
          else
          ax.set_ylabel(labels[i],fontsize="medium")

          end
      else
      end
    end
  end
 fig.align_ylabels()
 fig.align_xlabels()
  fig.autofmt_xdate()
end

# Corner plot for posterior distributions of 6 parameters, compared to true values
function corner(xs,true_xs,labels,nbins)
  x1,x2,x3,x4,x5,x6=xs[1],xs[2],xs[3],xs[4],xs[5],xs[6]
  truex1,truex2,truex3,truex4,truex5,truex6=true_xs
  label1,label2,label3,label4,label5,label6=labels

  fig=figure(figsize=(8,8))#,dpi=150)
  fig.subplots_adjust(hspace=0.25,wspace=0.25)
  ax1=fig.add_subplot(6,6,1)
  ax1.axvline(truex6,linestyle="--",color="black")
  ax1.hist(x6,bins=nbins,histtype="step",color="black")
  ax1.axvspan(quantile(vec(x6),0.1587),quantile(vec(x6),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
      #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false)
  ax2=fig.add_subplot(6,6,7,sharex=ax1)
  # ax2.hist2d(x6,x5,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=ax2.hist2d(x6,x5,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax2.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  ax2.set_ylabel(label5)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # ax2.minorticks_on()
  # ax2.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=true,labelbottom=false)
  ax3=fig.add_subplot(6,6,8)
  ax3.hist(x5,bins=nbins,histtype="step",color="black")
  ax3.axvline(truex5,linestyle="--",color="black")
  ax3.axvspan(quantile(vec(x5),0.1587),quantile(vec(x5),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # ax6.minorticks_on()
  # ax6.tick_params(which="both",direction="in",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false)
  ax4=fig.add_subplot(6,6,13,sharex=ax1)
  # ax3.hist2d(x6,x4,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=ax4.hist2d(x6,x4,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax4.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  ax4.set_ylabel(label4)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # ax9.minorticks_on()
  # ax9.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=true,labelbottom=false)
  ax5=fig.add_subplot(6,6,14,sharex=ax3)
  # ax4.hist2d(x5,x4,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=ax5.hist2d(x5,x4,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax5.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false) 
  # ax9.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=true,labelbottom=false)
  ax6=fig.add_subplot(6,6,15)
  ax6.hist(x4,bins=nbins,histtype="step",color="black")
  ax6.axvline(truex4,linestyle="--",color="black")
  ax6.axvspan(quantile(vec(x4),0.1587),quantile(vec(x4),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # ax5.minorticks_on()
  # ax5.tick_params(which="both",direction="in",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false) 
  ax7=fig.add_subplot(6,6,19,sharex=ax1)
  # ax6.hist2d(x6,x3,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=ax7.hist2d(x6,x3,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax7.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  ax7.set_ylabel(label3)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # ax6.minorticks_on()
  # ax6.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=false,labelbottom=false)
  ax8=fig.add_subplot(6,6,20,sharex=ax3)
  # ax7.hist2d(x5,x3,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=ax8.hist2d(x5,x3,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax8.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false) 
  # ax7.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=false,labelbottom=false)
  ax9=fig.add_subplot(6,6,21,sharex=ax6)
  # ax8.hist2d(x4,x3,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=ax9.hist2d(x4,x3,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax9.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # ax8.minorticks_on()
  # ax8.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=false,labelbottom=false)
  ax10=fig.add_subplot(6,6,22)
  ax10.hist(x3,bins=nbins,histtype="step",color="black")
  ax10.axvline(truex3,linestyle="--",color="black")
  ax10.axvspan(quantile(vec(x3),0.1587),quantile(vec(x3),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # ax9.minorticks_on()
  # ax9.tick_params(which="both",direction="in",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false) 
  ax11=fig.add_subplot(6,6,25,sharex=ax1)
  h,xedges,yedges=ax11.hist2d(x6,x2,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax11.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  ax11.set_ylabel(label2)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
      #     labelleft=false,labelbottom=false)
  ax12=fig.add_subplot(6,6,26,sharex=ax3)
  h,xedges,yedges=ax12.hist2d(x5,x2,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax12.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
      #     labelleft=false,labelbottom=false)
  ax13=fig.add_subplot(6,6,27,sharex=ax6)
  h,xedges,yedges=ax13.hist2d(x4,x2,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax13.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
      #     labelleft=false,labelbottom=false)
  ax14=fig.add_subplot(6,6,28,sharex=ax10)
  h,xedges,yedges=ax14.hist2d(x3,x2,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax14.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
      #     labelleft=false,labelbottom=false)
  ax15=fig.add_subplot(6,6,29)
  ax15.hist(x2,bins=nbins,histtype="step",color="black")
  ax15.axvline(truex2,linestyle="--",color="black")
  ax15.axvspan(quantile(vec(x2),0.1587),quantile(vec(x2),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
      #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false) 
  ax16=fig.add_subplot(6,6,31,sharex=ax1)
  h,xedges,yedges=ax16.hist2d(x6,x1,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax16.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  ax16.set_xlabel(label6)
  ax16.set_ylabel(label1)
  tick_params(left=false,bottom=false,labelbottom=false,labelleft=false)
  ax17=fig.add_subplot(6,6,32,sharex=ax3)
  h,xedges,yedges=ax17.hist2d(x5,x1,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax17.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  ax17.set_xlabel(label5)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # tick_params(axis="x",labelrotation=15)
      #     labelleft=false,labelbottom=true)
  ax18=fig.add_subplot(6,6,33,sharex=ax6)
  h,xedges,yedges=ax18.hist2d(x4,x1,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax18.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  ax18.set_xlabel(label4)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # tick_params(axis="x",labelrotation=15)
      #     labelleft=false,labelbottom=true)
  ax19=fig.add_subplot(6,6,34,sharex=ax10)
  h,xedges,yedges=ax19.hist2d(x3,x1,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax19.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  ax19.set_xlabel(label3)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # tick_params(axis="x",labelrotation=15)
      #     labelleft=false,labelbottom=true)
  ax20=fig.add_subplot(6,6,35,sharex=ax15)
  h,xedges,yedges=ax20.hist2d(x2,x1,bins=nbins,alpha=0.95)
  h_sort=sort(vec(h))
  sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
  sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
  c=ax20.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="white",levels=[h_sort[sig2[1]],h_sort[sig1[1]]])
  ax20.set_xlabel(label2)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # tick_params(axis="x",labelrotation=15)
      #     labelleft=false,labelbottom=true)
  ax21=fig.add_subplot(6,6,36)
  ax21.hist(x1,bins=nbins,histtype="step",color="black")
  ax21.axvline(truex1,linestyle="--",color="black")
  ax21.axvspan(quantile(vec(x1),0.1587),quantile(vec(x1),0.8413),color="limegreen",alpha=0.35)
  ax21.set_xlabel(label1)
  tick_params(left=false,labelleft=false,labelbottom=false,bottom=false)
  # tick_params(axis="x",labelrotation=15)
  fig.autofmt_xdate()
  fig.align_ylabels()
  fig.align_xlabels()
      #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=true,labeltop=false,labelleft=false,labelright=false)
  # tight_layout()
  # tick_params(left=false,labelleft=false,bottom=false,labelbottom=false)
 # return fig

end

function plot_traces()
  fig=figure(figsize=(8,6))
  fig,axs=plt.subplots(nrows=5,ncols=nplanet)
  fig.suptitle(string(model," traces; BIC=",BIC))
  for ax in axs
    for i=1:nparam
    for j=1:nwalkers 
    ax.plot(par_mcmc[j,iburn:end,i])
    end
    ax.set_ylabel(pname[i])
    end
  end
  title=string("IMAGES/trace/case",case_num,grid_type_nplanet,"-",sigma,"secs",nyear,"yrs.png")
  savefig(title)
end
# Create a corner plot for posterior distributions of planet parameters
function corner_plot(sigma,nyear,grid_type_nplanet,case_num,nbins,include_moon::Bool=false) 
  EM=true
  if case_num==1  #&& isfile(string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile=string("FITS/fromEMB/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
  elseif case_num==2 #&& isfile(string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    fitfile=string("FITS/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
  else
    return  println("MCMC file for case",case_num," with ",grid_type_nplanet," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  if include_moon
    EM=false
  end

  jldmc=jldopen(String(mcfile),"r")
  jldfit=jldopen(String(fitfile),"r")
  nwalkers,nsteps=jldmc["nwalkers"],jldmc["nsteps"]
  iburn,samples=jldmc["iburn"], jldmc["indepsamples"]
  par_mcmc=jldmc["par_mcmc"]; lprob_mcmc=jldmc["lprob_mcmc"]  ; param=jldmc["param"]
  pname=jldmc["pname"]
 
  # println("# of independent samples: ",samples)
  # println("Burn-in number: ",iburn," out of ",nsteps," steps")
  tt0,tt,ttmodel,sigtt=jldfit["tt0"],jldfit["tt"],jldfit["ttmodel"],jldfit["sigtt"]
  nplanet,ntrans=jldfit["nplanet"],jldfit["ntrans"]
  nt1,nt2=jldfit["ntrans"][1],jldfit["ntrans"][2]
  jmax=5
  weight=ones(nt1+nt2)./ sigtt.^2 
  nparam=length(par_mcmc[1,end,:])
  avg=zeros(nparam)
  minus1sig=zeros(nparam)
  plus1sig=zeros(nparam)
  for i=1:nparam
   avg[i],minus1sig[i],plus1sig[i]=quantile(vec(par_mcmc[:,iburn:end,i]),[0.5,0.1587,0.8413])
   # println(pname[i]," = ",avg[i]," + ",abs(plus1sig[i]-avg[i])," _ ",abs(avg[i]-minus1sig[i]))
  end
  # stable=false 
  # mutual_radius=mutual_Hill(avg[2],avg[1].*CGS.MSUN/CGS.MEARTtranspose(H),CGS.MSUN,avg[7],avg[6].*CGS.MSUN/CGS.MEARTH)
 # if mutual_radius <

  vals=jldmc["par_mcmc"][:,jldmc["iburn"]:end,:]#,sigdigits=6)
  BIC_from_mc,BIC_from_est=calc_BIC(jldmc["lprob_mcmc"][:,jldmc["iburn"]:jldmc["nsteps"]],jldfit["tt0"],jldfit["tt"],jldfit["sigtt"],jldfit["nplanet"],jldfit["ntrans"],vals,EM=EM)
  # println("BIC= ",BIC ,'\t'," reduced χ^2: ",reduced_chisq,'\t'," χ^2: ",chisq)

  sigsys=round((median(vec(par_mcmc[:,iburn:end,end]))).* 3600*24,sigdigits=3)
  sigsys_err=(std(vec(par_mcmc[:,iburn:end,end]))).* 3600*24
  sigtot=round(sqrt(sigsys^2 + sigma^2),sigdigits=4)

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
    L"$m_b / M_{⋆}$",L"$P_b$ [days]",L"$t_{0,b}$",L"$e_b cos(ω_b)$",L"$e_b sin(ω_b)$",
    L"$m_c / M_{⋆}$",L"$P_c$ [days]",L"$t_{0,c}$",L"$e_c cos(ω_c)$",L"$e_c sin(ω_c)$",
    L"$m_e / M_{⋆}$",L"$P_e$ [days]",L"$t_{0,e}$",L"$e_e cos(ω_e)$",L"$e_e sin(ω_e)$",
    L"$m_d / M_{⋆}$",L"$P_d$ [days]",L"$t_{0,d}$",L"$e_d cos(ω_d)$",L"$e_d sin(ω_d)$",
    L"$μ_f$",L"$P_f$ [days]",L"$t_{0,f}$",L"$e_f cos(ω_f)$",L"$e_f sin(ω_f)$",
    L"$t_{max} sin(ϕ_0)$",L"$t_{max} cos(ϕ_0)$",L"$Δϕ$ [rad]",L"$σ_{sys}^2$ [days]"]  

    # True values based on "PlanetaryBodyData.pdf" (source?)
    offset_p1=224.70
    m1=vec(par_mcmc[:,iburn:nsteps,1]).* CGS.MSUN/CGS.MEARTH
    ec1=vec(par_mcmc[:,iburn:nsteps,4])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
    es1=vec(par_mcmc[:,iburn:nsteps,5])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
    p1=vec(par_mcmc[:,iburn:nsteps,2])#.-offset_p1
    truem1=0.815
    trueec1=calc_evec1(0.00677323,131.53298)
    truees1=calc_evec2(0.00677323,131.53298)
    truep1=224.7007992#.-offset_p1
    # lim=0.00076,0.00081
    lim=minimum(p1),maximum(p1)#0.0064,0.00652
    label=L"Per$_1 - 224.7$ [days]"
    offset_p2=365.25
    m2=vec(par_mcmc[:,iburn:nsteps,6]).* CGS.MSUN/CGS.MEARTH;
    ec2=vec(par_mcmc[:,iburn:nsteps,9]);#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
    es2=vec(par_mcmc[:,iburn:nsteps,10]);#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
    p2=vec(par_mcmc[:,iburn:nsteps,7])#.-offset_p2;
    truem2=1
    trueec2=calc_evec1(0.01671022,102.94719)
    truees2=calc_evec2(0.01671022,102.94719)
    truep2=365.256355#-offset_p2 #365.256355

    values=[ 
    vec(par_mcmc[:,iburn:end,1]), 
    vec(par_mcmc[:,iburn:end,2]),
    vec(par_mcmc[:,iburn:end,3]), 
    vec(par_mcmc[:,iburn:end,6]), 
    vec(par_mcmc[:,iburn:end,7]), 
    vec(par_mcmc[:,iburn:end,8]),
    vec(par_mcmc[:,iburn:end,11]),
    vec(par_mcmc[:,iburn:end,12]),
     vec(par_mcmc[:,iburn:end,13])]
    # vec(par_mcmc[:,iburn:end,11]), vec(par_mcmc[:,iburn:end,12]),vec(par_mcmc[:,iburn:end,14]), vec(par_mcmc[:,iburn:end,15]),
    # vec(par_mcmc[:,iburn:end,16]), vec(par_mcmc[:,iburn:end,17]),vec(par_mcmc[:,iburn:end,19]), vec(par_mcmc[:,iburn:end,20])]
    labels=[parname[1];
    parname[2];
    parname[3];
    parname[6];
    parname[7];
    parname[8];
    parname[11];
    parname[12];
    parname[13]]
   # parname[11];parname[12];parname[14];parname[15];
   # parname[16];parname[17];parname[19];parname[20]]
    corner(values,labels,nbins;quantiles=[0.1587,0.8413])
    println("Plotting ",size(values)," values with ",size(labels),"labels.")
    # show()
    title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"-",sigma,"secs",nyear,"yrst0.png")
    # fig1=corner([m1,ec1,p1,m2,ec2,p2],[truem1,trueec1,truep1,truem2,trueec2,truep2],labels,nbins)
    # fig1.suptitle(string(model," Posteriors for Planet 1"))
    # fig1.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sigtot," sec",'\n',"BIC= ",BIC,'\n',L"$\chi^2 =$",chi))
    savefig(title)
    # clf()
end
    # if include_moon && grid_type_nplanet=="p3moon"
    #  tmax=sqrt.(vec(par_mcmc[:,iburn:nsteps,16]).^2 .+ vec(par_mcmc[:,iburn:nsteps,17]).^2)
    #  x1=vec(par_mcmc[:,iburn:nsteps,16])
    #  x2=vec(par_mcmc[:,iburn:nsteps,17])
    #  x3=vec(par_mcmc[:,iburn:nsteps,18])#.*57.2957795
    #  truetmax=calc_tmax(CGS.AU,CGS.AMOON*CGS.AU,CGS.MEARTtranspose(H),CGS.MMOON,365.256355) #0.0018
    #  truex2=0.01
    #  truex3=2.31586#.*57.2957795
    #  title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Moon2-",sigma,"secs",nyear,"yrs.png")
    #  fig5=corner(tmax,x3,truetmax,truex3,nbins)

    #  elseif include_moon 
    #    tmax=sqrt.(vec(par_mcmc[:,iburn:nsteps,21]).^2 .+ vec(par_mcmc[:,iburn:nsteps,22]).^2) .* 24*60
    #    x1=vec(par_mcmc[:,iburn:nsteps,21])
    #    x2=vec(par_mcmc[:,iburn:nsteps,22])
    #    x3=vec(par_mcmc[:,iburn:nsteps,23])#.*57.2957795
    #    truetmax=calc_tmax(CGS.AU,CGS.AMOON*CGS.AU,CGS.MEARTtranspose(H),CGS.MMOON,365.256355) .* 24*60#0.0018
    #    truex2=0.01
    #    truex3=2.31586#.*57.2957795
    #    title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Moon2-",sigma,"secs",nyear,"yrs.png")
    #    fig5=corner(tmax,x3,truetmax,truex3,nbins)
  #   title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"Moon2-",sigma,"secs",nyear,"yrs.png")
  #   corner(tmax,x3,truetmax,truex3,nbins)
  #   xlabel(L"$t_{max}$ [days]")
  #   ylabel(L"$\Delta \phi$ [rad]")
  #   tight_layout()
  #   savefig(title)
  #   clf()
  # end
# Basic corner plot for posterior distributions of 2 parameters
function corner(x,y,nbins)
  function scatter_hist(x, y, ax, ax_histx, ax_histy)
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=false)
    ax_histy.tick_params(axis="y", labelleft=false)
    # the scatter plot:
    # ax.scatter(x, y)
    h,xedges,yedges=ax.hist2d(x,y,bins=nbins,alpha=0.9)
    # the contour:
    # M=meshgrid(x,y)
    # xedges=h[2][2:end];yedges=h[3][2:end];
    h_sort=sort(vec(h))
    sig1=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.683))
    sig2=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.955))
    # sig3=findall((cumsum(h_sort)/sum(h_sort)).>= (1-0.997))
    # println(h_sort[sig1[1]])
    # ax.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),colors="black",levels=)
    ax.contour(xedges[1:end-1],yedges[1:end-1],transpose(h),levels=[h_sort[sig2[1]],h_sort[sig1[1]]],colors="white")
    # now determine nice limits by hand:
    # binwidth = 0.25
    # xymax = maximum(maximum(abs.(x)), maximum(abs.(y)))
    # lim = (int(xymax/binwidth) + 1) * binwidth
    # bins = range(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=nbins,histtype="step",density=true)
    ax_histy.hist(y, bins=nbins,histtype="step",density=true,orientation="horizontal")
  end
  fig=figure(figsize=(4,4))#,dpi=150)
  # gs = fig.add_gridspec(2, 2)
  # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
  # the size of the marginal axes and the main axes in both directions.
  gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
  # Create the Axes.
  ax = fig.add_subplot(gs[2, 1])
  ax_histx = fig.add_subplot(gs[1, 1], sharex=ax)
  ax_histy = fig.add_subplot(gs[2, 2], sharey=ax)
  # Draw the scatter plot and marginals.
  scatter_hist(x, y, ax, ax_histx, ax_histy)
    # return fig
end

function plot_acc(nyear,nplanet,grid_type_nplanet,set_ylim::Bool=false)
	# med90=ones(10).*NaN;med60=ones(10).*NaN
	med90,errors90=mc_vals(90,nyear,grid_type_nplanet,1)
	med60,errors60=mc_vals(60,nyear,grid_type_nplanet,1)
	med30,errors30=mc_vals(30,nyear,grid_type_nplanet,1)
	med10,errors10=mc_vals(10,nyear,grid_type_nplanet,1)
	labels=["planet b",
	"planet c",
	"planet e",
	"planet d"]
	colors=["salmon",
	"forestgreen",
	"orange",
	"firebrick"]
	true_x=[.815,
	1.012,
	0.1074,
	317.8]
	noise1=[10,10,10,10]
	noise3=[30,30,30,30]
	noise6=[60,60,60,60]
	noise9=[90,90,90,90]
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
	fig,ax=plt.subplots(figsize=(5,4))
	function plot_noise(ax)
	ax.set_title(string(model," for ",nyear," yr baseline"),fontsize="large")
	for iplanet=1:nplanet
	ax.errorbar(sqrt.(noise1[iplanet].^2 .+ (med10[end] .* 24 .*3600) .^2),med10[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors10[1,(iplanet-1)*5+1] errors10[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	ax.errorbar(sqrt.(noise3[iplanet].^2 .+ (med30[end] .* 24 .*3600) .^2),med30[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors30[1,(iplanet-1)*5+1] errors30[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	ax.errorbar(sqrt.(noise6[iplanet].^2 .+ (med60[end] .* 24 .*3600) .^2),med60[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors60[1,(iplanet-1)*5+1] errors60[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	###      ax.errorbar(sqrt.(noise9[iplanet].^2 .+ (med90[end] .* 24 .*3600) .^2),med90[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors90[1,(iplanet-1)*5+1] errors90[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	ax.axhline(true_x[iplanet],color="k",linestyle="--")
	end
	ax.set_ylabel(L"M$_p$ [$M_{\oplus}$]",fontsize="large")
	ax.set_xlabel(L"$\sqrt{\sigma_{obs}^2 +\sigma_{sys}^2}$ [secs]",fontsize="large")
	ax.minorticks_on()
	if set_ylim
	ax.set_ylim(-10,400)
	end
	ax.set_xlim(0,175)
	tight_layout()
	ax.text(45,317.8.*1.05,string(labels[4]),fontsize="medium")
	ax2=fig.add_axes([0.7,0.3,0.25,0.45])
	ax3=fig.add_axes([0.7,0.3])
	gs=fig.add_gridspec(2,1 ,left=0.6,right=0.9,bottom=0.3,top=0.8,wspace=0.05,hspace=0.05)
	ax1=fig.add_subplot(gs[2,1]);ax2=fig.add_subplot(gs[1,1],sharex=ax1)
	ax1.spines["top"].set_visible(false);ax2.spines["bottom"].set_visible(false)
	ax1.tick_params(labeltop=false);ax2.tick_params(labelbottom=false,bottom=false)
	d=0.02
	ax2.plot((-d, +d),(-d,+d),transform=ax2.transAxes,clip_on=false,linewidth=0.8,color="k")
	ax2.plot((1-d, 1+d),(-d,+d),transform=ax2.transAxes,clip_on=false,linewidth=0.8,color="k")
	ax1.plot((1-d, 1+d),(1-d,1+d),transform=ax1.transAxes,clip_on=false,linewidth=0.8,color="k")
	ax1.plot((-d, +d),(1-d,1+d),transform=ax1.transAxes,clip_on=false,linewidth=0.8,color="k")
	# ax2.text(10,true_x[1].*1.05,"b",fontsize="medium")
	# ax2.text(10,true_x[2].*1.025,"c",fontsize="medium")
	ax2.text(50,true_x[3].*1.5,string(labels[3]),fontsize="medium")
	for iplanet=1:3
	#for j=1:length(noise)
	ax2.errorbar(sqrt.(noise1[iplanet].^2 .+ (med10[end] .* 24 .*3600) .^2),med10[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors10[1,(iplanet-1)*5+1] errors10[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	ax2.errorbar(sqrt.(noise3[iplanet].^2 .+ (med30[end] .* 24 .*3600) .^2),med30[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors30[1,(iplanet-1)*5+1] errors30[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	ax2.errorbar(sqrt.(noise6[iplanet].^2 .+ (med60[end] .* 24 .*3600) .^2),med60[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors60[1,(iplanet-1)*5+1] errors60[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	### ax2.errorbar(sqrt.(noise9[iplanet].^2 .+ (med90[end] .* 24 .*3600) .^2),med90[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors90[1,(iplanet-1)*5+1] errors90[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	ax2.axhline(true_x[iplanet],color="k",linestyle="--")
	end


	#    iplanet=3
	#      ax1.errorbar(noise1[iplanet],med10[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors10[1,(iplanet-1)*5+1] errors10[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	#      ax1.errorbar(noise1[iplanet],med10[(iplanet-1)*5+1].*CGS.MSUN/CGS.MEARTH,yerr=[errors10[1,(iplanet-1)*5+1] errors10[2,(iplanet-1)*5+1]].*CGS.MSUN/CGS.MEARTH,fmt="o",color=colors[iplanet])
	ax2.minorticks_on();ax2.set_xlim(0,110);ax2.tick_params(which="both",direction="in");#ax2.set_ylim(0.01,1.1)
	end # plot_noise function
	plot_noise(ax)
	# savefig(string("IMAGES/discussion/",nplanet,"_case1",nyear,".png"),dpi=200)
end
function plot_param(grid_type_nplanet,nplanet::Real,case::Real=1,comparison_value::String="%err")
	sigmas=[10,30,60,90,120]
	nyears=[15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
	#   sigs=[];yrs=[]
	values_2d=zeros(length(sigmas),length(nyears)) .*NaN
	errors_2d=zeros(length(sigmas),length(nyears)) .*NaN
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
	fig=figure(figsize=(12,5))
	#gs=fig.add_gridspec(2,1)
	ax1=fig.add_subplot(1,nplanet,1)
	ax2=fig.add_subplot(1,nplanet,2)
	function make_plot(param_col,errors_2d,values_2d,ax)
	for (i,sigma) in enumerate(sigmas)
	for (j,nyear) in enumerate(nyears)
	mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
	if case==2
	mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
	end
	if isfile(mcfile)
	mc=jldopen(mcfile,"r")
	med=median( vec(mc["par_mcmc"][:,mc["iburn"]:end,param_col]))
	st_dev=std( vec(mc["par_mcmc"][:,mc["iburn"]:end,param_col]))
	av=mean( vec(mc["par_mcmc"][:,mc["iburn"]:end,param_col]))
	#                 low=quantile( vec(mc["par_mcmc"][:,mc["iburn"]:end,param_col]),0.1587)
	#                 high=quantile( vec(mc["par_mcmc"][:,mc["iburn"]:end,param_col]),0.8413)
	#                 worst=max(st_dev,med.-low,high.-med)
	#                 err_ratio=worst/med 
	#                push!(sigs,sigma);push!(yrs,nyear)
	if st_dev/av < 1
	errors_2d[i,j]=(st_dev/av)*100
	values_2d[i,j]=(abs(true_vals[param_col]-av)/st_dev)
	end
	#                 values_2d[i,j]=(med-true_vals[param_col])/true_vals[param_col]
	end
	end
	end
	parname=[
	L"$m_b / M_{⋆}$",L"$P_b$ [days]",L"$t_{0,b}$",L"$e_b cos(ω_b)$",L"$e_b sin(ω_b)$",
	L"$m_c / M_{⋆}$",L"$P_c$ [days]",L"$t_{0,c}$",L"$e_c cos(ω_c)$",L"$e_c sin(ω_c)$",
	#    L"$m_e / M_{⋆}$",L"$P_e$ [days]",L"$t_{0,e}$",L"$e_e cos(ω_e)$",L"$e_e sin(ω_e)$",
	L"$m_d / M_{⋆}$",L"$P_d$ [days]",L"$t_{0,d}$",L"$e_d cos(ω_d)$",L"$e_d sin(ω_d)$",
	L"$μ_f$",L"$P_f$ [days]",L"$t_{0,f}$",L"$e_f cos(ω_f)$",L"$e_f sin(ω_f)$",
	L"$t_{max} sin(ϕ_0)$",L"$t_{max} cos(ϕ_0)$",L"$Δϕ$ [rad]",L"$σ_{sys}^2$ [days]"]  

	#   fig,ax=plt.subplots(figsize=(6,5))
	#     h,xedges,yedges,im=ax.hist2d(yrs,sigs,weights=errors,bins=[length(yrs),length(sigs)],cmin=0.000001)
	if comparison_value=="%err"
	im=ax.imshow(errors_2d,aspect=3,cmap="plasma",origin="lower",vmin=0.00001,extent=[15,30,1,6])
	# im=ax.matshow(errors_2d,origin="lower",aspect=2) 
	ax.set_title(string("% Error in ",parname[param_col]," for ",model,),fontsize="large")
	elseif comparison_value=="true"
	im=ax.imshow(values_2d,aspect=1,cmap="plasma",origin="lower",vmin=0.00001)#,extent=[0.0,120,0.0,120])
	# ax.text(string(L"$\sigma$ "," of the Measurement of ",parname[param_col]," for ",model,),fontsize="large")
	end
	colorbar(im,ax=ax)
	ax.set_yticklabels(["0","10","30","60","90","120"])
	ax.set_xlabel(string("Baseline, ",L"$n_{year}$ [yrs]"),fontsize="large")
	end
	#     ylim=(0,120)
	ax1.set_ylabel(string(L"$\sigma_{obs}$ [secs]"),fontsize="large")
	make_plot(1,errors_2d,values_2d,ax1)
	values_2d=zeros(length(sigmas),length(nyears)) .*NaN
	errors_2d=zeros(length(sigmas),length(nyears)) .*NaN
	make_plot(6,errors_2d,values_2d,ax2)
	if nplanet==3
	ax3=fig.add_subplot(1,nplanet,3)
	values_2d=zeros(length(sigmas),length(nyears)) .*NaN
	errors_2d=zeros(length(sigmas),length(nyears)) .*NaN
	make_plot(11,errors_2d,values_2d,ax3)
	end

	# name=string("IMAGES/discussion/err","_",grid_type_nplanet,".jpg")  
	# savefig(name,dpi=200)
end


function comp_BIC(grid_type_nplanet,case=1,include_moon=false;grid_type_nplanet2::String="p2",grid_type_nplanet3::String="p3",grid_type_nplanet4::String="widep4")
	EM =true
  if include_moon
      EM=false
    end
    sigmas=[10,30,60,90]
    nyears=[15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    values_pn=ones(length(sigmas),length(nyears)) .* NaN
    values_p2=ones(length(sigmas),length(nyears)) .* NaN
    values_p3=ones(length(sigmas),length(nyears)) .* NaN
     values_p4=ones(length(sigmas),length(nyears)) .* NaN
	fig,axs=subplots(1,2,figsize=(10,5))
	images=[]
	#ax1=fig.add_subplot(121)
	#ax2=fig.add_subplot(122)
	function plot_BIC(grid_type_nplanet,values_pn,ax,fig)
	 if  grid_type_nplanet=="p2" 
	    model2=L"$\mathcal{H}_{PP}$"
	  elseif grid_type_nplanet=="p3" || grid_type_nplanet=="widep3"
	    model3=L"$\mathcal{H}_{PPP}$"
	  elseif grid_type_nplanet=="p4" || grid_type_nplanet=="widep4"
	    model4=L"$\mathcal{H}_{PPPP}$"
	  elseif grid_type_nplanet=="p3moon" || grid_type_nplanet=="widep3moon"
	    model=L"$\mathcal{H}_{PPsP}$"
	  elseif grid_type_nplanet=="p3moonp4" || grid_type_nplanet=="widep3moonp4"
	    model=L"$\mathcal{H}_{PPsPP}$"
	  end  

	    for (i,sigma) in enumerate(sigmas)
	        for (j,nyear) in enumerate(nyears)
			mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")

			mcfile2=string("MCMC/fromEMB/",grid_type_nplanet2,"_mcmc",sigma,"s",nyear,"yrs.jld2")
			mcfile3=string("MCMC/fromEMB/",grid_type_nplanet3,"_mcmc",sigma,"s",nyear,"yrs.jld2")

			fitfile=string("FITS/fromEMB/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
			fitfile4=string("FITS/fromEMB/",grid_type_nplanet4,"_fit",sigma,"s",nyear,"yrs.jld2")
			fitfile3=string("FITS/fromEMB/",grid_type_nplanet3,"_fit",sigma,"s",nyear,"yrs.jld2")
			fitfile2=string("FITS/fromEMB/",grid_type_nplanet2,"_fit",sigma,"s",nyear,"yrs.jld2")
	   if case==2
			mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
			fitfile=string("FITS/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
	    end
	   if isfile(fitfile) && isfile(mcfile) ##&& isfile(fitfile2) && isfile(fitfile)
			f=jldopen(String(fitfile),"r");f2=jldopen(String(fitfile2),"r");f3=jldopen(String(fitfile3),"r") 
			f=jldopen(String(fitfile),"r");f2=jldopen(String(fitfile2),"r");f3=jldopen(String(fitfile3),"r") 
			f4=jldopen(String(fitfile4),"r") 
			nplanet=f["nplanet"]
	# values_p2[i,j]=fit_BIC(f2["lprob_mcmc"],f2["tt0"],f2["tt"],f2["sigtt"],f2["nplanet"],f2["ntrans"],f2["best_p2"],EM=EM)[2]
	 values_p2[i,j]=marg_BIC(f2["lprob_best_p2"],f2["tt0"],f2["tt"],f2["sigtt"],f2["nplanet"],f2["ntrans"],f2["best_p2"],EM=EM)
	values_p3[i,j]=marg_BIC(f3["lprob_best_p3"],f3["tt0"],f3["tt"],f3["sigtt"],f3["nplanet"],f3["ntrans"],f3["best_p3"],EM=EM) 
		# values_p4[i,j]=fit_BIC(f4["lprob_best_p4"],f4["tt0"],f4["tt"],f4["sigtt"],f4["nplanet"],f4["ntrans"],f4["best_p4"],EM=EM)[2] 
	values_pn[i,j]=marg_BIC(f["lprob_best_p$nplanet"],f["tt0"],f["tt"],f["sigtt"],f["nplanet"],f["ntrans"],f["best_p$nplanet"],EM=EM)#,sigdigits=6)
	#end
	 # @show sigma,nyear,f["best_p5"][22]./365.25
	   end
	   end
	   end 
	   datasets=[[values_p3.-values_p2] [values_pn.-values_p3]]
	   models=[L"$\mathcal{H}_{PPP}$",L"$\mathcal{H}_{PPPP}$"]
	   norm=pltco.Normalize(vmin=minimum(minimum.(datasets)),vmax=maximum(maximum.(datasets)))
	   #@show datasets[:,:,:]
	   for i=1:2
	   ax[i].set_yticklabels(["0","10","30","60","90"])
	   ax[i].set_xlabel(string("Baseline, ",L"$n_{year}$ [yrs]"),fontsize="large")
	   ax[i].set_ylabel(string("Injected Noise, ",L"$\sigma_{obs}$ [s]"),fontsize="large")
	   # ax[i].set_title(string("Evidence for ",models[i]),fontsize="x-large")
	   # ax[i].imshow(values_pn,origin="lower",cmap="viridis",aspect=2,extent=[15,30,1,5])
	   push!(images,ax[i].imshow(datasets[i],origin="lower",cmap="viridis",aspect=2,extent=[15,30,1,5],interpolation="nearest",norm=norm))
	   end
	   #  im=ax[2].imshow(values_pn.-values_p3,aspect=2,cmap="viridis",origin="lower",extent=[15,30,1,6])
	#for i in 1:length(sigmas)
	#for j in 1:length(nyears)
	#if !isnan(values_pn[i,j])# 0.0
	#@show nyears[j],sigmas[i]
	#end
	#end
	#end
	return fig,images
	end # plot_BIC function
	fig,images=plot_BIC(grid_type_nplanet,values_pn,axs,fig)
	fig.subplots_adjust(hspace=0.05,wspace=0.22)
	fig.colorbar(images[1],ax=axs,fraction=0.1,label=L"$\Delta$BIC",orientation="horizontal")#,spacing="proportional")
	#vals_2=plot_BIC(grid_type_nplanet2,values_p2,ax2)
	#fig=figure(figsize=(6,5))
	#ax=fig.add_subplot(111)
	#ax2.set_title(L"$\Delta$ BIC")
	#im=ax2.imshow(abs.(vals_n).-abs.(vals_p2),aspect=2,cmap="viridis",origin="lower",extent=[15,30,1,6])
	# colorbar(im,fraction=0.05,orientation="horizontal",label=L"$\Delta$BIC")#,spacing="proportional")

	#   ax[i].text(nyears[j],sigmas[i],string(round(datasets[i],sigdigits=3)), ha="center",va="center",color="white")#,fontsize="x-small")
	      name=string("IMAGES/discussion/deltaBIC",case,"_",grid_type_nplanet,"-p3.jpg") 
	      #tight_layout()
	     #fig.show()
	    savefig(name,dpi=200)
	     return  fig
end

function hist_panel(values,nbins,labels)
	K=length(xs)
	K=length(values)
	fig,axs=plt.subplots(K,K+1,figsize=(8,3))
	for (i,x) in enumerate(values)
		for j=1:length(x)
		ax=axs[i,j]
		ax.hist(x[j],histtype="step",bins=nbins,label=labels[j],density=true)
		ax.set_xlabel(labels[j])
	   end
   end
	axs[1,1].set_ylabel("Hpp")
	axs[2,1].set_ylabel("Hppp")
	#     fig.legend()
	tight_layout()
	return fig
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
    L"$m_b / M_{⋆}$",L"$P_b$ [days]",L"$t_{0,b}$",L"$e_b cos(ω_b)$",L"$e_b sin(ω_b)$",
    L"$m_c / M_{⋆}$",L"$P_c$ [days]",L"$t_{0,c}$",L"$e_c cos(ω_c)$",L"$e_c sin(ω_c)$",
    L"$m_e / M_{⋆}$",L"$P_e$ [days]",L"$t_{0,e}$",L"$e_e cos(ω_e)$",L"$e_e sin(ω_e)$",
    L"$m_d / M_{⋆}$",L"$P_d$ [days]",L"$t_{0,d}$",L"$e_d cos(ω_d)$",L"$e_d sin(ω_d)$",
    L"$μ_f$",L"$P_f$ [days]",L"$t_{0,f}$",L"$e_f cos(ω_f)$",L"$e_f sin(ω_f)$"]
   
   mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
	mcfile2=string("MCMC/fromEMB/",grid_type_nplanet2,"_mcmc",sigma,"s",nyear,"yrs.jld2")
	mcfile3=string("MCMC/fromEMB/",grid_type_nplanet3,"_mcmc",sigma,"s",nyear,"yrs.jld2")
   if isfile(mcfile)
		jldmc=jldopen(mcfile,"r")

		nwalkers,nsteps=jldmc["nwalkers"],jldmc["nsteps"]
		iburn,samples=jldmc["iburn"], jldmc["indepsamples"]
		par_mcmc=jldmc["par_mcmc"]; lprob_mcmc=jldmc["lprob_mcmc"]  ; param=jldmc["param"]
		pname=jldmc["pname"]
		values=[ 
    vec(par_mcmc[:,iburn:end,1]), 
    vec(par_mcmc[:,iburn:end,2]),
    vec(par_mcmc[:,iburn:end,3])] 
    # vec(par_mcmc[:,iburn:end,5]), 
    # vec(par_mcmc[:,iburn:end,6]), 
    # vec(par_mcmc[:,iburn:end,7]),
    # vec(par_mcmc[:,iburn:end,9]),
    # vec(par_mcmc[:,iburn:end,10])]#,    vec(par_mcmc[:,iburn:end,11]), vec(par_mcmc[:,iburn:end,12]),vec(par_mcmc[:,iburn:end,14]), vec(par_mcmc[:,iburn:end,15]),
    #vec(par_mcmc[:,iburn:end,16]), vec(par_mcmc[:,iburn:end,17]),vec(par_mcmc[:,iburn:end,19]), vec(par_mcmc[:,iburn:end,20])]
    labels=[parname[1];
    parname[2];
    parname[3]];
    # parname[5];
    # parname[6];
    # parname[7];
    # parname[9];
    # parname[10]]#;   parname[11];parname[12];parname[14];parname[15];
#   parname[16];parname[17];parname[19];parname[20]]
    true_vals=[]
	hist_panel(values[1], true_vals, labels[1],nbins;model)
	end

end
