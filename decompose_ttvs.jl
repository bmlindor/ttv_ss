include("ttv_nplanet.jl")

function decompose_ttvs(nplanet, ntrans, params)
	jmax = 5
	pair_ttvs = zeros(nplanet,nplanet,maximum(ntrans))
	for i=1:nplanet-1, j=i+1:nplanet
		param = [params[(i-1)*5+1:i*5]; params[(j-1)*5+1:j*5]]
		ttv = ttv_nplanet(2, jmax, [ntrans[i];ntrans[j]], param)
		pair_ttvs[i,j,1:ntrans[i]] = ttv[1,1:ntrans[i]] #planet i wrt planet j
		pair_ttvs[j,i,1:ntrans[j]] = ttv[2,1:ntrans[j]] #planet j wrt planet i
	end
  return pair_ttvs
end