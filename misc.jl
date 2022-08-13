#include TTVFaster  
function find_coeffs(tt,period,sigtt)
    nt = length(tt)
    x = zeros(2,nt)
    x[1,1:nt] .= 1.0
    x[2,1] = 0.0 
    for i=2:nt
      x[2,i] = round((tt[i]-tt[1])/period) 
    end
    coeff,covcoeff = regress(x,tt,sigtt)
    return coeff,covcoeff
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


