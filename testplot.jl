using PyPlot,Statistics
rc("font",family="serif")

function corner_planet(x1,x2,x3,x4)
#     mu_1,P_1,t01,e1cosw1,e1sinw1,
# mu_2,P_2,t02,e2cosw2,e2sinw2,
# mu_3,P_3,t03,e3cosw3,e3sinw3,
# tmaxsinphi0,tmaxcosphi0,deltaphi = f["pbest_global"]
#     if np=1
        
        
	meanx=mean(xvalue);sigmax=std(xvalue)
	meany=mean(yvalue);sigmay=std(yvalue)

	fig=figure(figsize=(8,8))
	subplots_adjust(hspace=0.03,wspace=0.03)

	subplot(4,4,1)
	hist(x4,bins=nbins,histtype="step",density="true",color="black")
	minorticks_on()
	tick_params(which="major",direction="in",length=6,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")

    subplot(4,4,5)
    hist2D(x4,x3,bins=nbins,cmin=1)
    ylabel(L"$e \sin \varpi $")
    minorticks_on()
    tick_params(which="major",direction="in",top="true",right="true",length=6
    ,labelleft="true",labelbottom="false")
    tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="true",labelbottom="false")
    
    subplot(4,4,6)
    hist(x3,bins=nbins,histtype="step",density="true",color="black")
    minorticks_on()
	tick_params(which="major",direction="in",length=6,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
    
    subplot(4,4,9)
    hist2D(x4,x2,bins=nbins,cmin=1)
    ylabel(L"$e \cos \varpi $")
    minorticks_on()
    tick_params(which="major",direction="in",top="true",right="true",length=6
    ,labelleft="true",labelbottom="false")
    tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="true",labelbottom="false")
    
    subplot(4,4,10)
    hist2D(x3,x2,bins=nbins,cmin=1)
    minorticks_on()
    tick_params(which="major",direction="in",top="true",right="true",length=6
    ,labelleft="false",labelbottom="false")
    tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="false",labelbottom="false")
	
    subplot(4,4,11)
    hist(x2,bins=nbins,histtype="step",density="true",color="black")
    minorticks_on()
	tick_params(which="major",direction="in",length=6,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")
	tick_params(which="minor",direction="in",length=2,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="false",labeltop="false",labelleft="false",labelright="false")    
    
    subplot(4,4,13)
    hist2D(x4,x1,bins=nbins,cmin=1)
    xlabel("Per[days]")
    ylabel(L"Mass [$M_{\oplus}$]")
    minorticks_on()
    tick_params(which="major",direction="in",top="true",right="true",length=6
    ,labelleft="true")
    tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="true")
    
    subplot(4,4,14)
    hist2D(x3,x1,bins=nbins,cmin=1)
    xlabel(L"$e \sin \varpi $")
    minorticks_on()
    tick_params(which="major",direction="in",top="true",right="true",length=6
    ,labelleft="false")
    tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="false")
    
    subplot(4,4,15)
    hist2D(x2,x1,bins=nbins,cmin=1)
    xlabel(L"$e \cos \varpi $")
    minorticks_on()
    tick_params(which="major",direction="in",top="true",right="true",length=6
    ,labelleft="false",labelbottom="true")
    tick_params(which="minor",direction="in",top="true",right="true",length=2
    ,labelleft="false",labelbottom="true")
	
    subplot(4,4,16)
    hist(x1,bins=nbins,histtype="step",density="true",color="black")
    xlabel(L"Mass [$M_{\oplus}$]")
    minorticks_on()
	tick_params(which="major",direction="in",length=6,
	    left="false",right="false",top="true",bottom="true",
	    labelbottom="true",labeltop="false",labelleft="false",labelright="false")
	tick_params(which="minor",direction="in",length=2,
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

# 	subplot(223,sharex=ax2,sharey=ax3)
# 	ax1=gca()
# 	h = ax1.hexbin(xvalue,yvalue,bins=nbins,gridsize=nbins,mincnt=1)
# # 	plot(truex,truey,ms=8,mec="yellow",mfc="yellow",marker="s")
#     axvline(truex,linestyle="--",color="black")
#     axhline(truey,linestyle="--",color="black")
# 	ax1.axis([minimum(xvalue),maximum(xvalue),minimum(yvalue),maximum(yvalue)])
# 	ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
# 	ax1.tick_params(which="minor",direction="in",top="true",right="true",length=2)
	show()
end


