using TTVFaster,DataFrames,CSV,LsqFit
function chisquare(tt0,nplanet,ntrans,params,tt,sigtt,jmax,EM)
  chisq = 0.0  #check memory allocation >>>>>>>>>>>>
  # println(params,tt[1],sigtt[1])
  tt_model = ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM) 
  for j=1:length(tt)
    chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
  end
  return chisq
end

function second_peak_params(grid_file::String)
	data, header=readdlm(grid_file,',',header=true)
	pnts=length(data[:,1]) ; nparam=length(data[1,:])-1
	new_params=zeros(nparam)
	function scond_peak(data)
	for i=1:pnts
		if abs(data[i,end] - maximum(data[:,end]))<1 && data[i,end]!=maximum(data[:,end])
			#println("local maximum: ",data[i,end])
			return i
    end
  end
	end
	indx=scond_peak(data) ; new_params,lprob=data[indx,1:nparam],data[indx,end]
	#data=CSV.read(grid_file,DataFrame,header=true)
	#sorted=sort!(data,[:lprob])
	return new_params,lprob
end
