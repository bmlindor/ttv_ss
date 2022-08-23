#include TTVFaster  
function chisquare(tt0,nplanet,ntrans,params,tt,sigtt,jmax,EM)#,fixp3::Bool = false,p3_cur::Float64 = 0.0)
  chisq = 0.0  #check memory allocation >>>>>>>>>>>>
  # println(params,tt[1],sigtt[1])
  tt_model = ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM) 
  for j=1:length(tt)
    chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
  end
  # println(nplanet)
  return chisq
end

function second_peak(xgrid,lprob)
	second_peak=false
	s=sortperm(lprob)
	sorted_grid=xgrid[s];sorted_prob=lprob[s]
	#sorted_params=params[s]
	for i=1:length(xgrid)
		global_max = round(sorted_grid[end],sigdigits=2)
		local_max = round(sorted_grid[end-1],sigdigits=2)
		if abs(sorted_prob[end] - sorted_prob[end-1])< 5 && global_max!=local_max
			#println(local_max," =/= ",global_max)
			second_peak=true      
		end
	end
	return second_peak
end
