include("TTVFaster/ttv_nplanet.jl")
#  jmax = 5
#   # Call ttv_nplanet:
#   ttv = ttv_nplanet(nplanet, jmax, ntrans, params[1:5*nplanet])
#   # We measure transit times, not TTVs, so add back in the linear ephemeris:
#   t01 = params[3]
#   per1 = params[2]
#   ttv1 = collect(range(t01,stop = t01+per1*(n1-1),length = n1)) #this doesnt account for skipped transits
#   for i=1:n1
#     ttv1[i]+= ttv[1,i]
#   end
#   t02 = params[8]
#   per2 = params[7]
#   ttv2 = collect(range(t02,stop = t02+per2*(n2-1),length = n2))
# # Make plot of decomposed simulated transit timing variations in minutes for years observed
#   for i=1:n2
#     if EMB
#       ttv2[i] += ttv[2,i]
#     else
#       #tmax = param[end-2]; phi0 = param[end-1]; deltaphi = param[end]
#       ts = params[end-2] #tmax sinphi0
#       tc = params[end-1] #tmax cosphi0
#       deltaphi = params[end]
#       ttv2[i] += ttv[2,i] + ts*cos((i-1)*deltaphi) + tc*sin((i-1)*deltaphi)
#     end
#   end
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