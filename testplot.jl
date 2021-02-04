using PyPlot,Statistics
rc("font",family="serif")

function corner_planet(x1,x2,x3,x4,nbins,
    truex1,truex2,truex3,truex4)
#     mu_1,P_1,t01,e1cosw1,e1sinw1,
# mu_2,P_2,t02,e2cosw2,e2sinw2,
# mu_3,P_3,t03,e3cosw3,e3sinw3,
# tmaxsinphi0,tmaxcosphi0,deltaphi = f["pbest_global"]
#     if np=1
        
        
# 	meanx=mean(xvalue);sigmax=std(xvalue)
# 	meany=mean(yvalue);sigmay=std(yvalue)

	fig=figure(figsize=(8,8))
	subplots_adjust(hspace=0.03,wspace=0.03)
	subplot(4,4,1)
    ax1=gca()
	ax1.hist(x4,bins=nbins,histtype="step",density="true",color="black")
    ax1.axvline(truex4,linestyle="-",color="black")
	ax1.minorticks_on()
	ax1.tick_params(which="major",direction="in",length=5,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	ax1.tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")

    subplot(4,4,5)
    ax5=gca()
    ax5.hist2D(x4,x3,bins=nbins,cmin=1)
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
    
    subplot(4,4,9)
    ax9=gca()
    ax9.hist2D(x4,x2,bins=nbins,cmin=1)
    ylabel(L"$e \cos \varpi $")
    ax9.minorticks_on()
    ax9.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="true",labelbottom="false")
    ax9.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="true",labelbottom="false")
    
    subplot(4,4,10)
    ax10=gca()
    ax10.hist2D(x3,x2,bins=nbins,cmin=1)
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
    
    subplot(4,4,13)
    ax13=gca()
    ax13.hist2D(x4,x1,bins=nbins,cmin=1)
    xlabel("Per [days]")
    ylabel(L"Mass [$M_{\oplus}$]")
    ax13.minorticks_on()
    ax13.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="true")
    ax13.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="true")
    
    subplot(4,4,14)
    ax14=gca()
    ax14.hist2D(x3,x1,bins=nbins,cmin=1)
    xlabel(L"$e \sin \varpi $")
    ax14.minorticks_on()
    ax14.tick_params(which="major",direction="in",top="true",right="true",length=5
    ,labelleft="false")
    ax14.tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="false")
    
    subplot(4,4,15)
    ax15=gca()
    ax15.hist2D(x2,x1,bins=nbins,cmin=1)
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


