#include TTVFaster  
function if_second_peak(xgrid,lprob)
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


function global_fit