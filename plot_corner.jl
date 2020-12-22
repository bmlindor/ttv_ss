

function corner_plot(xvalue,yvalue,nbins,optx,opty,truex,truey)
	meanx=mean(xvalue);sigmax=std(xvalue)
	meany=mean(yvalue);sigmay=std(yvalue)

	fig=figure(figsize=(6,6))
	subplots_adjust(hspace=0.02,wspace=0.02)

	subplot(221)
	ax2=gca()
	ax2.hist(xvalue,bins=nbins,histtype="step",density="true",color="green")
	axvline(meanx-sigmax,color="grey",label=L"$1\sigma$ Confidence Limit")
	axvline(meanx+sigmax,color="grey") 
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
	ax3.hist(yvalue,bins=nbins,histtype="step",density="true",color="green",orientation="horizontal")
	axhline(meany-sigmay,color="grey")
	axhline(meany+sigmay,color="grey")
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
	ax1.hist2d(xvalue,yvalue,bins=nbins,cmin=1)
	plot(truex,truey,ms=8,mec="yellow",mfc="indigo",marker="o")
	ax1.axis([minimum(xvalue),maximum(xvalue),minimum(yvalue),maximum(yvalue)])
	ax1.tick_params(which="major",direction="in",top="true",right="true",length=6)
	ax1.tick_params(which="minor",direction="in",top="true",right="true",length=2)
	show()
end

