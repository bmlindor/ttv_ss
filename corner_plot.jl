using JLD2,PyPlot,Statistics,Distributions,LinearAlgebra
using3D() # Needed to create a 3D subplot
#include("CGS.jl")
include("MCMC.jl")
# rc("font",family="sans-serif")
# rc("lines",linewidth=1.5)
# Basic corner plot for posterior distributions of 2 parameters
function corner(x,y,nbins)
  function scatter_hist(x, y, ax, ax_histx, ax_histy)
    # no labels
    ax_histx.tick_params(axis="x", labelbottom="False")
    ax_histy.tick_params(axis="y", labelleft="False")
    # the scatter plot:
    # ax.scatter(x, y)
    h,xedges,yedges=ax.hist2d(x,y,bins=nbins,cmin=1,cmax=100000)
    # the contour:
    # M=meshgrid(x,y)

    # xedges=h[2][2:end];yedges=h[3][2:end];
    # ax.contour(xedges,yedges,h,levels=[10,30,50])
    # now determine nice limits by hand:
    # binwidth = 0.25
    # xymax = maximum(maximum(abs.(x)), maximum(abs.(y)))
    # lim = (int(xymax/binwidth) + 1) * binwidth
    # bins = range(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=nbins,histtype="step")
    ax_histy.hist(y, bins=nbins,histtype="step", orientation="horizontal")
  end
  fig=figure(figsize=(5,5))#,dpi=150)
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
    return fig
end
  # subplots_adjust(hspace=0.05,wspace=0.05)
  # ax2=subplot(221)
  # ax2.hist(x1,bins=nbins,histtype="step",color="black")
  # ax2.minorticks_on()
  # ax2.tick_params(which="both",direction="out",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  # ax3=subplot(224)
  # ax3.hist(x2,bins=nbins,histtype="step",color="black",orientation="horizontal")
  # ax3.minorticks_on()
  # ax3.tick_params(which="both",direction="out",
  #     left=true,right=true,top=false,bottom=false,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false)

        
# Corner plot for posterior distributions of 2 parameters, compared to true values
function corner(x1,x2,truex1,truex2,nbins)
	fig=figure(figsize=(4,4),dpi=150)
	# subplots_adjust(hspace=0.05,wspace=0.05)
  # ax2=fig.add_axes(221)
  # text(0.5,0.5,"test")
  ax1=fig.add_subplot(223)
  ax1.hist2d(x1,x2,bins=nbins,cmin=1,cmax=100000)
  ax1.axis([minimum(x1),maximum(x1),minimum(x2),maximum(x2)])
  xlabel(L"$t_{max}$ [min], TTV Amplitude on Planet 2")
  ylabel(L"$\Delta\phi$ [rad], $Per_2$ Phase Offset")
  ax1.tick_params(which="both",direction="out",top=true,right=true)

	ax2=fig.add_subplot(221)
	h1=ax2.hist(x1,bins=nbins,histtype="step",color="black")
  ax2.axvspan(quantile(x1,0.1587),quantile(x1,0.8413),color="limegreen",alpha=0.35)
	axvline(truex1,linestyle="--",color="black",label="True Value")
	ax2.minorticks_on()
	ax2.tick_params(which="both",direction="out",
	    left=false,right=false,top=false,bottom=true,
	    labelbottom=true,labeltop=false,labelleft=false,labelright=false)

	ax3=fig.add_subplot(224,sharey=ax1)
	h2=ax3.hist(x2,bins=nbins,histtype="step",color="black",orientation="horizontal")
  ax3.axhspan(quantile(x2,0.1587),quantile(x2,0.8413),color="limegreen",alpha=0.35)
	axhline(truex2,linestyle="--",color="black")
	ax3.minorticks_on()
	ax3.tick_params(which="both",direction="out",
	    left=false,right=true,top=false,bottom=false,
	    labelbottom=false,labeltop=false,labelleft=false,labelright=true)
  return fig
end

# Corner plot for posterior distributions of 6 parameters, compared to true values
function corner(xs,true_xs,labels,nbins)
  x1,x2,x3,x4,x5,x6=xs[1],xs[2],xs[3],xs[4],xs[5],xs[6]
  truex1,truex2,truex3,truex4,truex5,truex6=true_xs
  label1,label2,label3,label4,label5,label6=labels
  fig=figure(figsize=(8,8))#,dpi=150)
  fig.subplots_adjust(hspace=0.2,wspace=0.2)
  ax1=fig.add_subplot(6,6,1)
  ax1.axvline(truex6,linestyle="--",color="black")
  ax1.hist(x6,bins=nbins,histtype="step",color="black")
  ax1.axvspan(quantile(vec(x6),0.1587),quantile(vec(x6),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax2=fig.add_subplot(6,6,7,sharex=ax1)
  # ax2.hist2d(x6,x5,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x6,x5,nbins)
  c=ax2.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax2.set_ylabel(label5)
  tick_params(labelbottom=false)
  # ax2.minorticks_on()
  # ax2.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=true,labelbottom=false)

  ax3=fig.add_subplot(6,6,8)
  ax3.hist(x5,bins=nbins,histtype="step",color="black")
  ax3.axvline(truex5,linestyle="--",color="black")
  ax3.axvspan(quantile(vec(x5),0.1587),quantile(vec(x5),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false)
  # ax6.minorticks_on()
  # ax6.tick_params(which="both",direction="in",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax3=fig.add_subplot(6,6,13)
  # ax3.hist2d(x6,x4,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x6,x4,nbins)
  c=ax3.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax3.set_ylabel(label4)
  tick_params(labelbottom=false)
  # ax9.minorticks_on()
  # ax9.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=true,labelbottom=false)

  ax4=fig.add_subplot(6,6,14)
  # ax4.hist2d(x5,x4,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x5,x4,nbins)
  c=ax4.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
  # ax9.minorticks_on()
  # ax9.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=true,labelbottom=false)

  ax5=fig.add_subplot(6,6,15)
  ax5.hist(x4,bins=nbins,histtype="step",color="black")
  ax5.axvline(truex4,linestyle="--",color="black")
  ax5.axvspan(quantile(vec(x4),0.1587),quantile(vec(x4),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false)
  # ax5.minorticks_on()
  # ax5.tick_params(which="both",direction="in",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false) 

  ax6=fig.add_subplot(6,6,19,sharex=ax1)
  # ax6.hist2d(x6,x3,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x6,x3,nbins)
  c=ax6.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax6.set_ylabel(label3)
  tick_params(labelbottom=false)
  # ax6.minorticks_on()
  # ax6.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=false,labelbottom=false)

  ax7=fig.add_subplot(6,6,20)
  # ax7.hist2d(x5,x3,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x5,x3,nbins)
  c=ax7.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
  # ax7.minorticks_on()
  # ax7.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=false,labelbottom=false)

  ax8=fig.add_subplot(6,6,21)
  # ax8.hist2d(x4,x3,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x4,x3,nbins)
  c=ax8.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
  # ax8.minorticks_on()
  # ax8.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=false,labelbottom=false)

  ax9=fig.add_subplot(6,6,22)
  ax9.hist(x3,bins=nbins,histtype="step",color="black")
  ax9.axvline(truex3,linestyle="--",color="black")
  ax9.axvspan(quantile(vec(x3),0.1587),quantile(vec(x3),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false)
  # ax9.minorticks_on()
  # ax9.tick_params(which="both",direction="in",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false) 

  ax10=fig.add_subplot(6,6,25,sharex=ax1)
    h,xedges,yedges=np.histogram2d(x6,x2,nbins)
  c=ax10.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax10.set_ylabel(label2)
  tick_params(labelbottom=false)
      #     labelleft=false,labelbottom=false)

  ax11=fig.add_subplot(6,6,26)
    h,xedges,yedges=np.histogram2d(x5,x2,nbins)
  c=ax11.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=false)

  ax12=fig.add_subplot(6,6,27)
    h,xedges,yedges=np.histogram2d(x4,x2,nbins)
  c=ax12.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=false)

  ax13=fig.add_subplot(6,6,28)
    h,xedges,yedges=np.histogram2d(x3,x2,nbins)
  c=ax13.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=false)

  ax14=fig.add_subplot(6,6,29)
  ax14.hist(x2,bins=nbins,histtype="step",color="black")
  ax14.axvline(truex2,linestyle="--",color="black")
  ax14.axvspan(quantile(vec(x2),0.1587),quantile(vec(x2),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false) 

  ax15=fig.add_subplot(6,6,31,sharex=ax1)
    h,xedges,yedges=np.histogram2d(x6,x1,nbins)
  c=ax15.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax15.set_xlabel(label6)
  ax15.set_ylabel(label1)
      #     labelbottom=true,labelleft=true)

  ax16=fig.add_subplot(6,6,32)
    h,xedges,yedges=np.histogram2d(x5,x1,nbins)
  c=ax16.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax16.set_xlabel(label5)
  tick_params(left=false,labelleft=false)
      #     labelleft=false,labelbottom=true)

  ax17=fig.add_subplot(6,6,33)
    h,xedges,yedges=np.histogram2d(x4,x1,nbins)
  c=ax17.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax17.set_xlabel(label4)
  tick_params(left=false,labelleft=false)
      #     labelleft=false,labelbottom=true)

  ax18=fig.add_subplot(6,6,34)
    h,xedges,yedges=np.histogram2d(x3,x1,nbins)
  c=ax18.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax18.set_xlabel(label3)
  tick_params(left=false,labelleft=false)
      #     labelleft=false,labelbottom=true)

  ax19=fig.add_subplot(6,6,35)
    h,xedges,yedges=np.histogram2d(x2,x1,nbins)
  c=ax19.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax19.set_xlabel(label2)
  tick_params(left=false,labelleft=false)
      #     labelleft=false,labelbottom=true)

  ax20=fig.add_subplot(6,6,36)
  ax20.hist(x1,bins=nbins,histtype="step",color="black")
  ax20.axvline(truex1,linestyle="--",color="black")
  ax20.axvspan(quantile(vec(x1),0.1587),quantile(vec(x1),0.8413),color="limegreen",alpha=0.35)
  ax20.set_xlabel(label1)
  tick_params(left=false,labelleft=false)
      #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=true,labeltop=false,labelleft=false,labelright=false)
  # tight_layout()
  # tick_params(left=false,labelleft=false,bottom=false,labelbottom=false)
 # return fig
end
# Corner plot for posterior distributions of 8 parameters, compared to true values
function corner(xs,true_xs,labels,nbins,model::LaTeXString)
  x1,x2,x3,x4,x5,x6,x7,x8=xs[1],xs[2],xs[3],xs[4],xs[5],xs[6],xs[7],xs[8]
  truex1,truex2,truex3,truex4,truex5,truex6,truex7,truex8=true_xs
  label1,label2,label3,label4,label5,label6,label7,label8=labels
  fig=figure(figsize=(8,6),dpi=150)
  fig.subplots_adjust(hspace=0.2,wspace=0.2)
  # fig.suptitle(string("Posteriors for [30 yr span] with model: ",model))
  ax1=fig.add_subplot(8,8,1)
  ax1.axvline(truex8,linestyle="--",color="black")
  ax1.hist(x8,bins=nbins,histtype="step",color="black")
  ax1.axvspan(quantile(vec(x8),0.1587),quantile(vec(x8),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax2=fig.add_subplot(8,8,9,sharex=ax1)
  # ax2.hist2d(x6,x5,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x8,x7,nbins)
  c=ax2.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax2.set_ylabel(label7,fontsize="large")
  tick_params(labelbottom=false,labelsize="small")
  # ax2.minorticks_on()
  # ax2.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=true,labelbottom=false)

  ax3=fig.add_subplot(8,8,10)
  ax3.hist(x7,bins=nbins,histtype="step",color="black")
  ax3.axvline(truex7,linestyle="--",color="black")
  ax3.axvspan(quantile(vec(x7),0.1587),quantile(vec(x7),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelrotation=15,labelbottom=false)
  # ax6.minorticks_on()
  # ax6.tick_params(which="both",direction="in",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false)

  ax4=fig.add_subplot(8,8,17,sharex=ax1)
  # ax3.hist2d(x6,x4,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x8,x6,nbins)
  c=ax4.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax4.set_ylabel(label6,fontsize="large")
  tick_params(labelbottom=false,labelsize="small")
  # ax9.minorticks_on()
  # ax9.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=true,labelbottom=false)

  ax5=fig.add_subplot(8,8,18,sharex=ax3,sharey=ax4)
  # ax4.hist2d(x5,x4,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x7,x6,nbins)
  c=ax5.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
  # ax9.minorticks_on()
  # ax9.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=true,labelbottom=false)

  ax6=fig.add_subplot(8,8,19)
  ax6.hist(x6,bins=nbins,histtype="step",color="black")
  ax6.axvline(truex6,linestyle="--",color="black")
  ax6.axvspan(quantile(vec(x6),0.1587),quantile(vec(x6),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelrotation=15,labelbottom=false)
  # ax5.minorticks_on()
  # ax5.tick_params(which="both",direction="in",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false) 

  ax7=fig.add_subplot(8,8,25,sharex=ax1)
  # ax6.hist2d(x6,x3,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x8,x5,nbins)
  c=ax7.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax7.set_ylabel(label5,fontsize="large")
  tick_params(labelbottom=false,labelsize="small")
  # ax6.minorticks_on()
  # ax6.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=false,labelbottom=false)

  ax8=fig.add_subplot(8,8,26,sharex=ax3,sharey=ax7)
  # ax7.hist2d(x5,x3,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x7,x5,nbins)
  c=ax8.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
  # ax7.minorticks_on()
  # ax7.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=false,labelbottom=false)

  ax9=fig.add_subplot(8,8,27,sharex=ax6,sharey=ax7)
  # ax8.hist2d(x4,x3,bins=nbins,cmin=1,cmax=100000)
  h,xedges,yedges=np.histogram2d(x6,x5,nbins)
  c=ax9.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
  # ax8.minorticks_on()
  # ax8.tick_params(which="both",direction="out",top=true,right=false,
  #     labelleft=false,labelbottom=false)

  ax10=fig.add_subplot(8,8,28)
  ax10.hist(x5,bins=nbins,histtype="step",color="black")
  ax10.axvline(truex5,linestyle="--",color="black")
  ax10.axvspan(quantile(vec(x5),0.1587),quantile(vec(x5),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelrotation=15,labelbottom=false)
  # ax9.minorticks_on()
  # ax9.tick_params(which="both",direction="in",
  #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false) 

  ax11=fig.add_subplot(8,8,33,sharex=ax1)
    h,xedges,yedges=np.histogram2d(x8,x4,nbins)
  c=ax11.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax11.set_ylabel(label4,fontsize="large")
  tick_params(labelbottom=false,labelsize="small")
      #     labelleft=false,labelbottom=false)

  ax12=fig.add_subplot(8,8,34,sharex=ax3,sharey=ax11)
    h,xedges,yedges=np.histogram2d(x7,x4,nbins)
  c=ax12.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=false)

  ax13=fig.add_subplot(8,8,35,sharex=ax6,sharey=ax11)
    h,xedges,yedges=np.histogram2d(x6,x4,nbins)
  c=ax13.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=false)

  ax14=fig.add_subplot(8,8,36,sharex=ax10,sharey=ax11)
    h,xedges,yedges=np.histogram2d(x5,x4,nbins)
  c=ax14.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=false)

  ax16=fig.add_subplot(8,8,37)
  ax16.hist(x4,bins=nbins,histtype="step",color="black")
  ax16.axvline(truex4,linestyle="--",color="black")
  ax16.axvspan(quantile(vec(x4),0.1587),quantile(vec(x4),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelrotation=15,labelbottom=false)
      #     left=false,right=false,top=true,bottom=true,
  #     labelbottom=false,labeltop=false,labelleft=false,labelright=false) 

  ax17=fig.add_subplot(8,8,41,sharex=ax1)
  h,xedges,yedges=np.histogram2d(x8,x3,nbins)
  c=ax17.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax17.set_ylabel(label3,fontsize="large")
  tick_params(labelbottom=false,labelsize="small")

      #     labelbottom=true,labelleft=true)

  ax18=fig.add_subplot(8,8,42,sharex=ax3,sharey=ax17)
  h,xedges,yedges=np.histogram2d(x7,x3,nbins)
  c=ax18.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  # ax18.set_xlabel(label7)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=true)

  ax19=fig.add_subplot(8,8,43,sharex=ax6,sharey=ax17)
    h,xedges,yedges=np.histogram2d(x6,x3,nbins)
  c=ax19.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=true)

  ax20=fig.add_subplot(8,8,44,sharex=ax10,sharey=ax17)
    h,xedges,yedges=np.histogram2d(x5,x3,nbins)
  c=ax20.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=true)

  ax21=fig.add_subplot(8,8,45,sharex=ax16,sharey=ax17)
    h,xedges,yedges=np.histogram2d(x4,x3,nbins)
  c=ax21.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=true)

  ax22=fig.add_subplot(8,8,46)
  ax22.hist(x3,bins=nbins,histtype="step",color="black")
  ax22.axvline(truex3,linestyle="--",color="black")
  ax22.axvspan(quantile(vec(x3),0.1587),quantile(vec(x3),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelrotation=15,labelbottom=false)

  ax23=fig.add_subplot(8,8,49,sharex=ax1)
    h,xedges,yedges=np.histogram2d(x8,x2,nbins)
  ax23.set_ylabel(label2,fontsize="large")
  c=ax23.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
      #     labelbottom=true,labelleft=true)
    tick_params(labelbottom=false,labelsize="small")

  ax24=fig.add_subplot(8,8,50,sharex=ax3,sharey=ax23)
    h,xedges,yedges=np.histogram2d(x7,x2,nbins)
  c=ax24.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=true)

  ax25=fig.add_subplot(8,8,51,sharex=ax6,sharey=ax23)
    h,xedges,yedges=np.histogram2d(x6,x2,nbins)
  c=ax25.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=true)

  ax26=fig.add_subplot(8,8,52,sharex=ax10,sharey=ax23)
    h,xedges,yedges=np.histogram2d(x5,x2,nbins)
  c=ax26.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=true)

  ax27=fig.add_subplot(8,8,53,sharex=ax16,sharey=ax23)
    h,xedges,yedges=np.histogram2d(x4,x2,nbins)
  c=ax27.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=true)
  
  ax28=fig.add_subplot(8,8,54,sharex=ax22,sharey=ax23)
    h,xedges,yedges=np.histogram2d(x3,x2,nbins)
  c=ax28.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  tick_params(left=false,labelleft=false,labelbottom=false)
      #     labelleft=false,labelbottom=true)

  ax29=fig.add_subplot(8,8,55)
  ax29.hist(x2,bins=nbins,histtype="step",color="black")
  ax29.axvline(truex2,linestyle="--",color="black")
  ax29.axvspan(quantile(vec(x2),0.1587),quantile(vec(x2),0.8413),color="limegreen",alpha=0.35)
  tick_params(left=false,labelleft=false,labelrotation=15,labelbottom=false)

  ax30=fig.add_subplot(8,8,57,sharex=ax1)
    h,xedges,yedges=np.histogram2d(x8,x1,nbins)
  c=ax30.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax30.set_xlabel(label8,fontsize="large")
  ax30.set_ylabel(label1,fontsize="large")
  tick_params(axis="x",labelrotation=20,labelsize="small")
  tick_params(axis="y",labelsize="small")
      #     labelbottom=true,labelleft=true)

  ax31=fig.add_subplot(8,8,58,sharex=ax3,sharey=ax30)
    h,xedges,yedges=np.histogram2d(x7,x1,nbins)
  c=ax31.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax31.set_xlabel(label7,fontsize="large")
  tick_params(left=false,labelleft=false,labelrotation=20,labelsize="small")
  # ax31.set_ylabel(label1)
      #     labelbottom=true,labelleft=true)

  ax32=fig.add_subplot(8,8,59,sharex=ax6,sharey=ax30)
    h,xedges,yedges=np.histogram2d(x6,x1,nbins)
  c=ax32.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax32.set_xlabel(label6,fontsize="large")
  tick_params(left=false,labelleft=false,labelrotation=20,labelsize="small")
      #     labelleft=false,labelbottom=true)

  ax33=fig.add_subplot(8,8,60,sharex=ax10,sharey=ax30)
    h,xedges,yedges=np.histogram2d(x5,x1,nbins)
  c=ax33.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax33.set_xlabel(label5,fontsize="large")
  tick_params(left=false,labelleft=false,labelrotation=20,labelsize="small")

      #     labelleft=false,labelbottom=true)
  ax34=fig.add_subplot(8,8,61,sharex=ax16,sharey=ax30)
    h,xedges,yedges=np.histogram2d(x4,x1,nbins)
  c=ax34.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax34.set_xlabel(label4,fontsize="large")
  tick_params(left=false,labelleft=false,labelrotation=20,labelsize="small")
      #     labelleft=false,labelbottom=true)
  
  ax35=fig.add_subplot(8,8,62,sharex=ax22,sharey=ax30)
  h,xedges,yedges=np.histogram2d(x3,x1,nbins)
  c=ax35.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax35.set_xlabel(label3,fontsize="large")
  tick_params(left=false,labelleft=false,labelrotation=20,labelsize="small")
      #     labelleft=false,labelbottom=true)

  ax36=fig.add_subplot(8,8,63,sharex=ax29,sharey=ax30)
    h,xedges,yedges=np.histogram2d(x2,x1,nbins)
  c=ax36.contour(xedges[1:end-1],yedges[1:end-1],h,cmap="viridis",normalize=true)
  ax36.set_xlabel(label2,fontsize="large")
  tick_params(left=false,labelleft=false,labelrotation=20,labelsize="small")
      #     labelleft=false,labelbottom=true)

  ax37=fig.add_subplot(8,8,64)
  ax37.hist(x1,bins=nbins,histtype="step",color="black")
  ax37.axvline(truex1,linestyle="--",color="black")
  ax37.axvspan(quantile(vec(x1),0.1587),quantile(vec(x1),0.8413),color="limegreen",alpha=0.35)
  ax37.set_xlabel(label1,fontsize="large")
  tick_params(left=false,labelleft=false,labelrotation=20,labelsize="small")
  fig.align_ylabels()
  fig.align_xlabels()
  # tight_layout()

 # return fig
end

# Create a corner plot for significant posterior distributions of planet parameters
function corner_hist(sigma,nyear,grid_type_nplanet,case_num,nbins,include_moon::Bool=false) 
  if case_num==1  && isfile(string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  elseif case_num==2 && isfile(string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
    mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  else
    return  println("MCMC file for case: ",case_num," with ",grid_type_nplanet," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
  end
  m=jldopen(String(mcfile),"r")
  par_mcmc= m["par_mcmc"]
  lprob_mcmc=m["lprob_mcmc"]
  nwalkers=m["nwalkers"]
  nsteps=m["nsteps"]
  iburn=m["iburn"]
  indepsamples=m["indepsamples"]
  # True values based on "PlanetaryBodyData.pdf" (source?)
  offset=224.70
  m1=vec(par_mcmc[:,iburn:nsteps,1]).* CGS.MSUN/CGS.MEARTH
  ec1=vec(par_mcmc[:,iburn:nsteps,4])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
  es1=vec(par_mcmc[:,iburn:nsteps,5])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
  p1=vec(par_mcmc[:,iburn:nsteps,2]).-offset
  e1=sqrt.(vec(par_mcmc[:,iburn:nsteps,4]).^2 .+ vec(par_mcmc[:,iburn:nsteps,5]).^2)
  truem1=0.815
  trueec1=calc_evec1(0.00677323,131.53298)
  truees1=calc_evec2(0.00677323,131.53298)
  truep1=224.7007992.-offset
  truee1=0.00677323
  lim=0.00076,0.00081
  offset=365.25
  m2=vec(par_mcmc[:,iburn:nsteps,6]).* CGS.MSUN/CGS.MEARTH
  ec2=vec(par_mcmc[:,iburn:nsteps,9])
  es2=vec(par_mcmc[:,iburn:nsteps,10])#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
  p2=vec(par_mcmc[:,iburn:nsteps,7]).-offset
  e2=sqrt.(vec(par_mcmc[:,iburn:nsteps,9]).^2 .+ vec(par_mcmc[:,iburn:nsteps,10]).^2)
  truem2=1
  trueec2=calc_evec1(0.01671022,102.94719)
  truees2=calc_evec2(0.01671022,102.94719)
  truep2=365.2564-offset #365.256355
  truee2=0.01671022
  lim=0.0064,0.00652
  corner(m1,m2,truem1,truem2,nbins)
  ylabel(L"Mass of Earth [$M_{Earth}$]")
  xlabel(L"Mass of Venus [$M_{Earth}$]")
  tight_layout()
  title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"masses",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  corner(e1,e2,truee1,truee2,nbins)
  ylabel("Eccentricity of Earth")
  xlabel("Eccentricity of Venus")
  tight_layout()
  title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"eccs",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  corner(ec1,ec2,trueec1,trueec2,nbins)
  ylabel(L"$e \cos \varpi$ for Earth")
  xlabel(L"$e \cos \varpi$ for Venus")
  tight_layout()
  title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"ecos",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  corner(es1,es2,truees1,truees2,nbins)
  ylabel(L"$e \sin \varpi$ for Earth")
  xlabel(L"$e \sin \varpi$ for Venus")
  tight_layout()
  title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"esin",sigma,"secs",nyear,"yrs.png")
  savefig(title)
  clf()
  if grid_type_nplanet=="p4"
    m3=vec(par_mcmc[:,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH
    ec3=vec(par_mcmc[:,iburn:nsteps,14])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    es3=vec(par_mcmc[:,iburn:nsteps,15])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    p3=vec(par_mcmc[:,iburn:nsteps,12])
    e3=sqrt.(vec(par_mcmc[:,iburn:nsteps,14]).^2 .+ vec(par_mcmc[:,iburn:nsteps,15]).^2)
    truem3=0.1074
    trueec3=calc_evec1(0.09341233,336.04084)
    truees3=calc_evec2(0.09341233,336.04084)
    truep3=686.9795859
    truee3=0.09341233
    corner(m3,e3,truem3,truee3,nbins)
    xlabel(L"Mass of Mars [$M_{Earth}$]")
    ylabel("Eccentricity of Mars")
    tight_layout()
    title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"Vmecc",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
    m4=vec(par_mcmc[:,iburn:nsteps,16]).* CGS.MSUN/CGS.MEARTH
    ec4=vec(par_mcmc[:,iburn:nsteps,19])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    es4=vec(par_mcmc[:,iburn:nsteps,20])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    p4=vec(par_mcmc[:,iburn:nsteps,17])
    e4=sqrt.(vec(par_mcmc[:,iburn:nsteps,19]).^2 .+ vec(par_mcmc[:,iburn:nsteps,20]).^2)
    truem4=317.8
    trueec4=calc_evec1(0.04839266,14.75385)
    truees4=calc_evec2(0.04839266,14.75385)
    truep4=4332.82012875
    truee4=0.04839266
    corner(m4,e4,truem4,truee4,nbins)
    xlabel(L"Mass of Jupiter [$M_{Earth}$]")
    ylabel("Eccentricity of Jupiter")
    tight_layout()
    title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"Jmecc",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
  else 
    m3=vec(par_mcmc[:,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH
    ec3=vec(par_mcmc[:,iburn:nsteps,14])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    es3=vec(par_mcmc[:,iburn:nsteps,15])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
    p3=vec(par_mcmc[:,iburn:nsteps,12])
    e3=sqrt.(vec(par_mcmc[:,iburn:nsteps,14]).^2 .+ vec(par_mcmc[:,iburn:nsteps,15]).^2)
    truem3=317.8
    trueec3=calc_evec1(0.04839266,14.75385)
    truees3=calc_evec2(0.04839266,14.75385)
    truep3=4332.82012875
    truee3=0.04839266
    corner(m3,e3,truem3,truee3,nbins)
    xlabel(L"Mass of Jupiter [$M_{Earth}$]")
    ylabel("Eccentricity of Jupiter")
    tight_layout()
    title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"mecc",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
  end
  if include_moon
    tmax=sqrt.(vec(par_mcmc[:,iburn:nsteps,16]).^2 .+ vec(par_mcmc[:,iburn:nsteps,17]).^2)
    x1=vec(par_mcmc[:,iburn:nsteps,16])
    x2=vec(par_mcmc[:,iburn:nsteps,17])
    x3=vec(par_mcmc[:,iburn:nsteps,18])#.*57.2957795
    truetmax=calc_tmax(CGS.AU,CGS.AMOON*CGS.AU,CGS.MEARTH,CGS.MMOON,365.256355) #0.0018
    truex2=0.01
    truex3=2.31586#.*57.2957795
    title=string("IMAGES/discussion/case",case_num,grid_type_nplanet,"Moon2-",sigma,"secs",nyear,"yrs.png")
    corner(tmax,x3,truetmax,truex3,nbins)
    xlabel(L"$t_{max}$ [days]")
    ylabel(L"$\Delta \phi$ [rad]")
    tight_layout()
    savefig(title)
    clf()
  end
  # show()
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
  # mutual_radius=mutual_Hill(avg[2],avg[1].*CGS.MSUN/CGS.MEARTH,CGS.MSUN,avg[7],avg[6].*CGS.MSUN/CGS.MEARTH)
 # if mutual_radius <
 # Calculate Bayesian Inference Criterion (BIC) 
  function calc_BIC(prob)
    function calc_chisq(par_mcmc)
    chisq = 0.0  
    tt_model = TTVFaster.ttv_wrapper(tt0,nplanet,ntrans,par_mcmc[1:nparam-1],jmax,EM) 
      for j=1:length(tt)
        chisq += (tt[j]-tt_model[j])^2 / (sigtt[j]^2 + par_mcmc[end]^2)
      end
    return chisq
    end
    chisq=calc_chisq(avg)
    N=length(tt0) ; k=nparam
    #println("[N_obs]= ",N," [no. of model params]= ",k)
    reduced_chisq=chisq/(N-k)
    BIC_chi(chisq,k,N)=chisq + k*log(N)
    BIC=-2*log(prob) + k*log(N)
    return reduced_chisq, BIC
  end
  prob=quantile(exp.(lprob_mcmc[iburn:nsteps]),0.5);prob_max = maximum(exp.(lprob_mcmc[iburn:nsteps]))
  #println(" median Prob: ",prob,"      maximum Prob: ",prob_max)
  #chi2_avg = chi_mcmc(tt0,nplanet,ntrans,mean_posteriors,tt,sigtt,jmax,EM)
  chi,BIC=round.(calc_BIC(prob_max),sigdigits=4)
  println("BIC= ",BIC ,'\t'," reduced χ^2: ",chi)

  sigsys=round((median(vec(par_mcmc[:,iburn:end,end]))).* 3600*24,sigdigits=3)
  sigsys_err=(std(vec(par_mcmc[:,iburn:end,end]))).* 3600*24
  sigtot=round(sqrt(sigsys^2 + sigma^2),sigdigits=4)

  #BIC,sigsys,sigtot,chi=mc_vals(sigma,nyear,grid_type_nplanet,case_num,include_moon)
  # println("Loaded",mcfile,".")
  # figure(figsize=(5,5))
  # plot(vec(par_mcmc[:,iburn:end,end]),vec(lprob_mcmc[:,iburn:end]))
  # xlabel(L"$σ_{sys}$")
  # ylabel(L"$\log(P)$")
  # @show

  # Find Percentage of walkers where difference between median and quantile value is >100
  bad_walk=[]
  for i in 1:nwalkers
    for j in 1:nparam
      walker_med,walker_quant=quantile!(par_mcmc[i,jldmc["iburn"]+1:end,j],[0.5,0.9])
      walk_start=par_mcmc[i,jldmc["iburn"]+1,j] 
      walk_end = par_mcmc[i,jldmc["iburn"]+1,j]
      ratio = walk_end/walk_start
      walker_prob=median(lprob_mcmc[i,jldmc["iburn"]+1:end])
      if abs(walk_end-walk_start)/walk_start > 0.1
        #abs(walker_med-walker_end)>30
        # println(i," ",walker_prob[i])
        append!(bad_walk,i)
      end
    end
  # If prob for a given chain is low, reject it

  #     # If systematic uncertainty > injected uncertainty, reject
  #   # if median(par_mcmc[i,jldmc["iburn"]:end,end]).*3600*24 >= sigma
  #   #   # println("Reject results?")
  #   #   append!(bad_walk,i)
  #   # end
  # end
  # println("Bad walkers: ",bad_walk)

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
  parname=[L"$μ_1$ [$M_{⋆}$]",L"$P_1$ [days]",L"$t_{0,1}$",L"$e_1 cos(ω_1)$",L"$e_1 sin(ω_1)$",
    L"$μ_2$ [$M_{⋆}$]",L"$P_2$ [days]",L"$t_{0,2}$",L"$e_2 cos(ω_2)$",L"$e_2 sin(ω_2)$",
    L"$μ_3$ [$M_{⋆}$]",L"$P_3$ [days]",L"$t_{0,3}$",L"$e_3 cos(ω_3)$",L"$e_3 sin(ω_3)$",
    L"$μ_4$ [$M_{⋆}$]",L"$P_4$ [days]",L"$t_{0,4}$",L"$e_4 cos(ω_4)$",L"$e_4 sin(ω_4)$",
    L"$μ_5$ [$M_{⋆}$]",L"$P_5$ [days]",L"$t_{0,5}$",L"$e_4 cos(ω_5)$",L"$e_5 sin(ω_5)$",
    L"$t_{max} sin(ϕ_0)$",L"$t_{max} cos(ϕ_0)$",L"$Δϕ$ [rad]",L"$σ_{sys}^2$ [days]"]
  function plot_traces()
    # fig=figure(figsize=(8,6))
    # fig,axs=plt.subplots(nrows=5,ncols=nplanet)
    # fig.suptitle(string(model," traces; BIC=",BIC))
    # for ax in axs
    #   for i=1:nparam
    #   for j=1:nwalkers 
    #   ax.plot(par_mcmc[j,iburn:end,i])
    #   end
    #   ax.set_ylabel(pname[i])
    #   end
    # end
    # title=string("IMAGES/trace/",grid_type_nplanet,"-",sigma,"secs",nyear,"yrs.png")
    # savefig(title)
    figure(figsize=(8,6))
    for i=1:5
    ax1=subplot(3,2,i)
    for j=1:nwalkers 
    ax1.plot(par_mcmc[j,iburn:nsteps,i])
    end
    ax1.set_ylabel(parname[i])
    end
    tight_layout()
    title=string("IMAGES/trace/",grid_type_nplanet,"Venus-",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()
    figure(figsize=(8,6))
    for i=1:5
      ax2=subplot(3,2,i)
      # for i=1:nparam
      # ax=subplot(nplanet,5,i)
      for j=1:nwalkers 
      ax2.plot(par_mcmc[j,iburn:nsteps,i+5])
      # ax.plot(par_mcmc[j,iburn:end,i])
      end
      ax2.set_ylabel(parname[i+5])
    end
    tight_layout()
    title=string("IMAGES/trace/",grid_type_nplanet,"Earth-",sigma,"secs",nyear,"yrs.png")
    savefig(title)
    clf()

    if nplanet==5
      figure(figsize=(8,6))
      for i=1:5
      ax3=subplot(3,2,i)
      for j=1:nwalkers 
      ax3.plot(par_mcmc[j,iburn:nsteps,i+20])
      end
      ax3.set_ylabel(parname[i+20])
      end
      tight_layout()
      title=string("IMAGES/trace/",grid_type_nplanet,"Saturn-",sigma,"secs",nyear,"yrs.png")
      savefig(title)
      clf()
      figure(figsize=(8,6))
      for i=1:5
      ax3=subplot(3,2,i)
      for j=1:nwalkers 
      ax3.plot(par_mcmc[j,iburn:nsteps,i+10])
      end
      ax3.set_ylabel(parname[i+10])
      end
      tight_layout()
      title=string("IMAGES/trace/",grid_type_nplanet,"Mars-",sigma,"secs",nyear,"yrs.png")
      savefig(title)
      #    println("Hit return to continue")
      #    read(STDIN,Char)
      clf()
      figure(figsize=(8,6))
      for i=1:5
      ax3=subplot(3,2,i)
      for j=1:nwalkers 
      ax3.plot(par_mcmc[j,iburn:nsteps,i+15])
      end
      ax3.set_ylabel(parname[i+15])
      end
      ax4=subplot(3,2,6)
      for j=1:nwalkers
        ax4.plot(par_mcmc[j,iburn:nsteps,end])
        ax4.set_ylabel(parname[end])
      end
      tight_layout()
      title=string("IMAGES/trace/",grid_type_nplanet,"Jupiter-",sigma,"secs",nyear,"yrs.png")
      savefig(title)
      clf()
    elseif nplanet==4
      figure(figsize=(8,6))
      for i=1:5
        ax3=subplot(3,2,i)
        for j=1:nwalkers 
        ax3.plot(par_mcmc[j,iburn:nsteps,i+10])
        end
        ax3.set_ylabel(parname[i+10])
      end
      tight_layout()
      title=string("IMAGES/trace/",grid_type_nplanet,"Mars-",sigma,"secs",nyear,"yrs.png")
      savefig(title)
      clf()
      figure(figsize=(8,6))
      for i=1:5
        ax3=subplot(3,2,i)
        for j=1:nwalkers 
        ax3.plot(par_mcmc[j,iburn:nsteps,i+15])
        end
        ax3.set_ylabel(parname[i+15])
      end
      ax4=subplot(3,2,6)
      for j=1:nwalkers
        ax4.plot(par_mcmc[j,iburn:nsteps,end])
        ax4.set_ylabel(parname[end])
      end
      tight_layout()
      title=string("IMAGES/trace/",grid_type_nplanet,"Jupiter-",sigma,"secs",nyear,"yrs.png")
      savefig(title)
      clf()
    else
      figure(figsize=(8,6))
      for i=1:5
      ax3=subplot(3,2,i)
      for j=1:nwalkers 
      ax3.plot(par_mcmc[j,iburn:nsteps,i+10])
      end
      ax3.set_ylabel(parname[i+10])
      end
      ax4=subplot(3,2,6)
      for j=1:nwalkers
        ax4.plot(par_mcmc[j,iburn:nsteps,end])
        ax4.set_ylabel(parname[end])
      end
      tight_layout()
      title=string("IMAGES/trace/",grid_type_nplanet,"Jupiter-",sigma,"secs",nyear,"yrs.png")
      savefig(title)
      clf()
    end
    if include_moon
      figure(figsize=(7,5))
      for i=1:3
      ax5=subplot(2,2,i)
      for j=1:nwalkers 
      ax5.plot(par_mcmc[j,1:nsteps,i+15])
      end
      ax5.set_ylabel(parname[i+20])
      end
      # subplot(2,2,4)
      # for j=1:nwalkers
      # plot(lprob_mcmc[j,iburn:nsteps])  
      # ylabel(L"$logProb$")
      # end
      tight_layout()
      title=string("IMAGES/trace/",grid_type_nplanet,"Moon-",sigma,"secs",nyear,"yrs.png")
      savefig(title)
      clf()
    end
  end
  #   for i=1:nparam
  #     for j=1:nwalkers 
  #     # ax.plot(par_mcmc[j,iburn:end,i])
  #     end
  #     # ax.set_ylabel(pname[i])
  #     # tight_layout()
  #     title=string("IMAGES/trace/",grid_type_nplanet,"-",sigma,"secs",nyear,"yrs.png")
  #     savefig(title)
  #     #    println("Hit return to continue")
  #     #    read(STDIN,Char)
  #     clf()
  #   end
  # end
  function plot_dist()
    # True values based on "PlanetaryBodyData.pdf" (source?)
    offset_p1=224.70
    m1=vec(par_mcmc[:,iburn:nsteps,1]).* CGS.MSUN/CGS.MEARTH
    ec1=vec(par_mcmc[:,iburn:nsteps,4])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
    es1=vec(par_mcmc[:,iburn:nsteps,5])#.*sqrt.(vec(par_mcmc[11:20,iburn:nsteps,4]).^2 .+ vec(par_mcmc[11:20,iburn:nsteps,5]).^2)
    p1=vec(par_mcmc[:,iburn:nsteps,2]).-offset_p1
    truem1=0.815
    trueec1=calc_evec1(0.00677323,131.53298)
    truees1=calc_evec2(0.00677323,131.53298)
    truep1=224.7007992.-offset_p1
    # lim=0.00076,0.00081
    lim=minimum(p1),maximum(p1)#0.0064,0.00652
    label=L"Per$_1 - 224.7$ [days]"
    title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Venus-",sigma,"secs",nyear,"yrs.png")
    fig1=corner(m1,ec1,es1,p1,truem1,trueec1,truees1,truep1,nbins,lim,label)
    fig1.suptitle(string(model," Posteriors for Planet 1"))
    fig1.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sigtot," sec",'\n',"BIC= ",BIC,'\n',L"$\chi^2 =$",chi))
    savefig(title)
    clf()

    offset_p2=365.25
    m2=vec(par_mcmc[:,iburn:nsteps,6]).* CGS.MSUN/CGS.MEARTH;
    ec2=vec(par_mcmc[:,iburn:nsteps,9]);#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
    es2=vec(par_mcmc[:,iburn:nsteps,10]);#.*sqrt.(vec(par_mcmc[:,1:nsteps,9]).^2 .+ vec(par_mcmc[:,1:nsteps,10]).^2)
    p2=vec(par_mcmc[:,iburn:nsteps,7]).-offset_p2;
    truem2=1
    trueec2=calc_evec1(0.01671022,102.94719)
    truees2=calc_evec2(0.01671022,102.94719)
    truep2=365.256355-offset_p2 #365.256355
    lim=minimum(p2),maximum(p2)#0.0064,0.00652
    label=L"Per$_2 - 365.25$ [days]"
    title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Earth-",sigma,"secs",nyear,"yrs.png")
    fig2=corner(m2,ec2,es2,p2,truem2,trueec2,truees2,truep2,nbins,lim,label)
    fig2.suptitle(string(model," Posteriors for Planet 2"))
    fig2.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sigtot," sec",'\n',"BIC= ",BIC,'\n',L"$\chi^2 =$",chi))
    savefig(title)
    clf()
    if grid_type_nplanet=="p4" || grid_type_nplanet=="p3moonp4" || grid_type_nplanet=="widep3moonp4" || grid_type_nplanet=="widep4"
      m3=vec(par_mcmc[:,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH
      ec3=vec(par_mcmc[:,iburn:nsteps,14])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
      es3=vec(par_mcmc[:,iburn:nsteps,15])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
      p3=vec(par_mcmc[:,iburn:nsteps,12])
      truem3=0.1074
      trueec3=calc_evec1(0.09341233,336.04084)
      truees3=calc_evec2(0.09341233,336.04084)
      truep3=686.9795859
      lim=minimum(p3),maximum(p3)
      label=L"Per$_4$ [days]"
      title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Mars-",sigma,"secs",nyear,"yrs.png")
      fig3=corner(m3,ec3,es3,p3,truem3,trueec3,truees3,truep3,nbins,lim,label)
      fig3.suptitle(string(model," Posteriors for Planet 4"))
      fig3.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sigtot," sec",'\n',"BIC= ",BIC,'\n',L"$\chi^2 =$",chi))
      savefig(title)
      clf()
      #  println("Mars= ",minimum(m3)," ",minimum(p3)," ",minimum(ec3)," ",minimum(es3))
      #  println("Mars= ",maximum(m3)," ",maximum(p3)," ",maximum(ec3)," ",maximum(es3))
      m4=vec(par_mcmc[:,iburn:nsteps,16]).* CGS.MSUN/CGS.MEARTH
      ec4=vec(par_mcmc[:,iburn:nsteps,19])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
      es4=vec(par_mcmc[:,iburn:nsteps,20])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
      p4=vec(par_mcmc[:,iburn:nsteps,17])
      truem4=317.8
      trueec4=calc_evec1(0.04839266,14.75385)
      truees4=calc_evec2(0.04839266,14.75385)
      truep4=4332.82012875
      lim=minimum(p4),maximum(p4)
      label=L"Per$_3$ [days]"
      title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Jupiter-",sigma,"secs",nyear,"yrs.png")
       # println("Jupiter= ",minimum(m4)," ",minimum(p4)," ",minimum(ec4)," ",minimum(es4))
       # println("Jupiter= ",maximum(m4)," ",maximum(p4)," ",maximum(ec4)," ",maximum(es4))
      fig4=corner(m4,ec4,es4,p4,truem4,trueec4,truees4,truep4,nbins,lim,label)
      fig4.suptitle(string(model," Posteriors for Planet 3"))
      fig4.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sigtot," sec",'\n',"BIC= ",BIC,'\n',L"$\chi^2 =$",chi))
      savefig(title)
      clf()
    elseif grid_type_nplanet=="p3" || grid_type_nplanet=="widep3" || grid_type_nplanet=="p3moon"
      m3=vec(par_mcmc[:,iburn:nsteps,11]).* CGS.MSUN/CGS.MEARTH
      ec3=vec(par_mcmc[:,iburn:nsteps,14])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
      es3=vec(par_mcmc[:,iburn:nsteps,15])#.*sqrt.(vec(par_mcmc[:,1:nsteps,14]).^2 .+ vec(par_mcmc[:,1:nsteps,15]).^2)
      p3=vec(par_mcmc[:,iburn:nsteps,12])
      truem3=317.8
      trueec3=calc_evec1(0.04839266,14.75385)
      truees3=calc_evec2(0.04839266,14.75385)
      truep3=4332.82012875
      lim=minimum(p3),maximum(p3)
      label=L"Per$_3$ [days]"
      title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Jupiter-",sigma,"secs",nyear,"yrs.png")
       # println("Jupiter= ",minimum(m3)," ",minimum(p3)," ",minimum(ec3)," ",minimum(es3))
       # println("Jupiter= ",maximum(m3)," ",maximum(p3)," ",maximum(ec3)," ",maximum(es3))
      fig3=corner(m3,ec3,es3,p3,truem3,trueec3,truees3,truep3,nbins,lim,label)
      fig3.suptitle(string(model," Posteriors for Planet 3"))
      fig3.text(0.36,0.8,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sigtot," sec",'\n',"BIC= ",BIC,'\n',L"$\chi^2 =$",chi))
      savefig(title)
      clf()
    end
    if include_moon && grid_type_nplanet=="p3moon"
     tmax=sqrt.(vec(par_mcmc[:,iburn:nsteps,16]).^2 .+ vec(par_mcmc[:,iburn:nsteps,17]).^2)
     x1=vec(par_mcmc[:,iburn:nsteps,16])
     x2=vec(par_mcmc[:,iburn:nsteps,17])
     x3=vec(par_mcmc[:,iburn:nsteps,18])#.*57.2957795
     truetmax=calc_tmax(CGS.AU,CGS.AMOON*CGS.AU,CGS.MEARTH,CGS.MMOON,365.256355) #0.0018
     truex2=0.01
     truex3=2.31586#.*57.2957795
     title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Moon2-",sigma,"secs",nyear,"yrs.png")
     fig5=corner(tmax,x3,truetmax,truex3,nbins)
     fig5.suptitle(string(model," Posteriors for Satellite"))
     fig5.text(0.575,0.725,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sigtot," sec",'\n',"BIC= ",BIC,'\n',L"$\chi^2 =$",chi))
      # corner(x1,x2,x3,nbins)
     savefig(title)
     clf()
     elseif include_moon 
       tmax=sqrt.(vec(par_mcmc[:,iburn:nsteps,21]).^2 .+ vec(par_mcmc[:,iburn:nsteps,22]).^2) .* 24*60
       x1=vec(par_mcmc[:,iburn:nsteps,21])
       x2=vec(par_mcmc[:,iburn:nsteps,22])
       x3=vec(par_mcmc[:,iburn:nsteps,23])#.*57.2957795
       truetmax=calc_tmax(CGS.AU,CGS.AMOON*CGS.AU,CGS.MEARTH,CGS.MMOON,365.256355) .* 24*60#0.0018
       truex2=0.01
       truex3=2.31586#.*57.2957795
       title=string("IMAGES/corner/case",case_num,grid_type_nplanet,"Moon2-",sigma,"secs",nyear,"yrs.png")
       fig5=corner(tmax,x3,truetmax,truex3,nbins)
       fig5.suptitle(string(model," Posteriors for Satellite"))
       fig5.text(0.575,0.725,string(L"$\sigma_{sys}=$",sigsys," sec",'\n',L"$\sigma_{tot}=$",sigtot," sec",'\n',"BIC= ",BIC,'\n',L"$\chi^2 =$",chi))
        #corner(x1,x2,x3,nbins)
       savefig(title)
       clf()
     end
  end
  plot_traces()
  # plot_dist()
end
