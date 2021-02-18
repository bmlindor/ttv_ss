using PyPlot,Statistics
rc("font",family="serif")

function corner_plot(xvalue,yvalue,nbins,optx,opty,truex,truey)
	meanx=mean(xvalue);sigmax=std(xvalue)
	meany=mean(yvalue);sigmay=std(yvalue)

	fig=figure(figsize=(6,6))
	subplots_adjust(hspace=0.02,wspace=0.02)

	subplot(221)
	ax2=gca()
	ax2.hist(xvalue,bins=nbins,histtype="step",density="true",color="black")
	axvline(meanx-sigmax,color="grey",alpha=0.5,label=L"$1\sigma$ Limit")
	axvline(meanx+sigmax,color="grey",alpha=0.5) 
	axvline(truex,linestyle="--",color="black",label="True Value")
	# plot(truex,truey,color="blue",marker="o",label="True Value")
	
	ax2.minorticks_on()
	ax2.tick_params(which="major",direction="in",length=6,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax2.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax2.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.0)

	subplot(224)
	ax3=gca()
	ax3.hist(yvalue,bins=nbins,histtype="step",density="true",color="black",orientation="horizontal")
	axhline(meany-sigmay,color="grey",alpha=0.5)
	axhline(meany+sigmay,color="grey",alpha=0.5)
	axhline(truey,linestyle="--",color="black")
	# axhline(truey,color="blue")
	ax3.minorticks_on()
	ax3.tick_params(which="major",direction="in",length=6,
	    left="true",right="true",top="false",bottom="false",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax3.tick_params(which="minor",direction="in",length=2,
	    left="true",right="true",top="false",bottom="false",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")

	subplot(223,sharex=ax2,sharey=ax3)
	ax1=gca()
	# h = ax1.hexbin(xvalue,yvalue,bins=nbins,gridsize=nbins,mincnt=1)
	ax1.hist2d(xvalue,yvalue,bins=nbins,cmin=1)
# 	plot(truex,truey,ms=8,mec="yellow",mfc="yellow",marker="s")
    axvline(truex,linestyle="--",color="black")
    axhline(truey,linestyle="--",color="black")
	ax1.axis([minimum(xvalue),maximum(xvalue),minimum(yvalue),maximum(yvalue)])
	ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
	ax1.tick_params(which="minor",direction="in",top="true",right="true",length=2)
	show()
end

function corner_moon(x1,x2,x3,nbins,
    truex1,truex2,truex3,pname)
    
    fig=figure(figsize=(8,8))
	subplots_adjust(hspace=0.05,wspace=0.05)
	subplot(3,3,1)
    ax1=gca()
    ax1.axvline(truex3,linestyle="-",color="black")
	ax1.hist(x3,bins=nbins,histtype="step",density="true",color="black")
	ax1.minorticks_on()
	ax1.tick_params(which="major",direction="in",length=5,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax1.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")

    subplot(3,3,4,sharex=ax1)
    ax4=gca()
    ax4.hist2d(x3,x2,bins=nbins,cmin=1)
    ylabel(L"$t_{max} \cos{\phi_0}$")
    ax4.minorticks_on()
    ax4.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="true",labelbottom="false")
    ax4.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="true",labelbottom="false")
    
    subplot(3,3,5)
    ax5=gca()
    ax5.hist(x2,bins=nbins,histtype="step",density="true",color="black")
#     ax5.axvline(truex2,linestyle="-",color="black")
    ax5.minorticks_on()
	ax5.tick_params(which="major",direction="in",length=5,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax5.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
    
    subplot(3,3,7,sharex=ax1)
    ax7=gca()
    ax7.hist2d(x3,x1,bins=nbins,cmin=1)
    xlabel(L"$\Delta \phi [deg]$")
    ylabel(L"$t_{max} \sin{\phi_0}$")
    ax7.minorticks_on()
    ax7.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="true",labelbottom="true")
    ax7.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="true",labelbottom="true")
    
    subplot(3,3,8,sharex=ax5)
    ax8=gca()
    ax8.hist2d(x2,x1,bins=nbins,cmin=1)
    xlabel(L"$t_{max} \cos{\phi_0}$")
    ax8.minorticks_on()
    ax8.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="false",labelbottom="true")
    ax8.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="false",labelbottom="true")
    
    subplot(3,3,9)
    ax9=gca()
    ax9.hist(x1,bins=nbins,histtype="step",density="true",color="black")
    xlabel(L"$t_{max} \sin{\phi_0}$")
#     ax9.axvline(truex1,linestyle="-",color="black")
#     ylabel(L"$e \cos \varpi $")
    ax9.minorticks_on()
    ax9.tick_params(which="major",direction="in",top="true",right="false",length=5
    ,labelleft="false",labelbottom="true")
    ax9.tick_params(which="minor",direction="in",top="true",right="false",length=2
    ,labelleft="false",labelbottom="true")
    show()
end


function corner_planet(x1,x2,x3,x4,nbins,
    truex1,truex2,truex3,truex4,pname)
#     mu_1,P_1,t01,e1cosw1,e1sinw1,
# mu_2,P_2,t02,e2cosw2,e2sinw2,
# mu_3,P_3,t03,e3cosw3,e3sinw3,
# tmaxsinphi0,tmaxcosphi0,deltaphi = f["pbest_global"]
#     if np=1
    if string(pname) == "venus"
        offset = 224.7007
    elseif string(pname) == "earth"
        offset = 365.2554
    else 
        offset = 0.0
    end
#     println(truex4-offset)
# 	meanx=mean(xvalue);sigmax=std(xvalue)
# 	meany=mean(yvalue);sigmay=std(yvalue)

	fig=figure(figsize=(8,8))
	subplots_adjust(hspace=0.05,wspace=0.05)
	subplot(4,4,1)
    ax1=gca()
    ax1.axvline(truex4-offset,linestyle="-",color="black")
	ax1.hist(x4.-offset,bins=nbins,histtype="step",density="true",color="black")
	ax1.minorticks_on()
	ax1.tick_params(which="major",direction="in",length=5,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax1.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")

    subplot(4,4,5,sharex=ax1)
    ax5=gca()
    ax5.hist2d(x4.-offset,x3,bins=nbins,cmin=1)
    ylabel(L"$e \sin \varpi $")
    ax5.minorticks_on()
    ax5.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="true",labelbottom="false")
    ax5.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="true",labelbottom="false")
    
    subplot(4,4,6)
    ax6=gca()
    ax6.hist(x3,bins=nbins,histtype="step",density="true",color="black")
    ax6.axvline(truex3,linestyle="-",color="black")
    ax6.minorticks_on()
	ax6.tick_params(which="major",direction="in",length=5,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax6.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
    
    subplot(4,4,9,sharex=ax1)
    ax9=gca()
    ax9.hist2d(x4.-offset,x2,bins=nbins,cmin=1)
    ylabel(L"$e \cos \varpi $")
    ax9.minorticks_on()
    ax9.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="true",labelbottom="false")
    ax9.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="true",labelbottom="false")
    
    subplot(4,4,10,sharex=ax6)
    ax10=gca()
    ax10.hist2d(x3,x2,bins=nbins,cmin=1)
    ax10.minorticks_on()
    ax10.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="false",labelbottom="false")
    ax10.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="false",labelbottom="false")
	
    subplot(4,4,11)
    ax11=gca()
    ax11.hist(x2,bins=nbins,histtype="step",density="true",color="black")
    ax11.axvline(truex2,linestyle="-",color="black")
    ax11.minorticks_on()
	ax11.tick_params(which="major",direction="in",length=5,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax11.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")    
    
    subplot(4,4,13,sharex=ax1)
    ax13=gca()
    ax13.hist2d(x4.-offset,x1,bins=nbins,cmin=1)
    xlabel(L"$Per - P_{offset}$ [days]")
    ylabel(L"Mass [$M_{\oplus}$]")
    ax13.minorticks_on()
    ax13.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelbottom="true",labelleft="true")
    ax13.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelbottom="true",labelleft="true")
    
    subplot(4,4,14,sharex=ax6)
    ax14=gca()
    ax14.hist2d(x3,x1,bins=nbins,cmin=1)
    xlabel(L"$e \sin \varpi $")
    ax14.minorticks_on()
    ax14.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="false")
    ax14.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="false")
    
    subplot(4,4,15,sharex=ax11)
    ax15=gca()
    ax15.hist2d(x2,x1,bins=nbins,cmin=1)
    xlabel(L"$e \cos \varpi $")
    ax15.minorticks_on()
    ax15.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="false",labelbottom="true")
    ax15.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="false",labelbottom="true")
	
    subplot(4,4,16)
    ax16=gca()
    ax16.hist(x1,bins=nbins,histtype="step",density="true",color="black")
    ax16.axvline(truex1,linestyle="-",color="black")
    xlabel(L"Mass [$M_{\oplus}$]")
    ax16.minorticks_on()
	ax16.tick_params(which="major",direction="in",length=5,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="true",labeltop="false",labelleft="false",labelright="false")
	ax16.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="true",labeltop="false",labelleft="false",labelright="false")
# 	ax3.hist(yvalue,bins=nbins,histtype="step",density="true",color="black",orientation="horizontal")
# # 	axhline(meany-sigmay,color="grey",alpha=0.5)
# # 	axhline(meany+sigmay,color="grey",alpha=0.5)
# # 	axhline(truey,linestyle="--",color="black")
# 	# axhline(truey,color="blue")
# 	ax3.minorticks_on()
# 	ax3.tick_params(which="major",direction="in",length=6,
# 	    left="true",right="true",top="false",bottom="false",
# 	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
# 	ax3.tick_params(which="minor",direction="in",length=2,
# 	    left="true",right="true",top="false",bottom="false",
# 	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	show()
end
