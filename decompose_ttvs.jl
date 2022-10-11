using TTVFaster

function decompose_ttvs(nplanet,ntrans,params)
	jmax = 5
	pair_ttvs = zeros(nplanet,nplanet,maximum(ntrans))
	for i=1:nplanet-1,j=i+1:nplanet
		param = [params[(i-1)*5+1:i*5]; params[(j-1)*5+1:j*5]]
		ttv = TTVFaster.ttv_nplanet(2,jmax,[ntrans[i];ntrans[j]],param)
		#ttv = ttv_nplanet(nplanet,jmax,[ntrans[i];ntrans[j]],[params[(i-1)*5+1:i*5]; params[(j-1)*5+1:j*5]])
		pair_ttvs[i,j,1:ntrans[i]] = ttv[1,1:ntrans[i]] #planet i wrt planet j
		pair_ttvs[j,i,1:ntrans[j]] = ttv[2,1:ntrans[j]] #planet j wrt planet i
	end
  return pair_ttvs
end

function moon_ttvs(ntrans,params)
	ttvs = zeros(ntrans[2])
   # We measure transit times,not TTVs,so add back in the linear ephemeris:
	for i=1:ntrans[2]
		  #tmax = param[end-2]; phi0 = param[end-1]; deltaphi = param[end]
		  ts = params[end-2] #tmax sinphi0
		  tc = params[end-1] #tmax cosphi0
		  deltaphi = params[end]
		  ttvs[i] = ts*cos((i-1)*deltaphi) + tc*sin((i-1)*deltaphi)
	end 
	return ttvs
end
