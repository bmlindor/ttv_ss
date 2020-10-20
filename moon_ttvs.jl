include("TTVFaster/ttv_nplanet.jl")

function moon_ttvs(nplanet, ntrans, params)
	jmax = 5
	moon_ttvs = zeros(maximum(ntrans))
	ttv1 = collect(range(0,stop=ntrans[1]-1,length=ntrans[1]))
	ttv2 = collect(range(0,stop=ntrans[2]-1,length=ntrans[2]))
	ttv = ttv_nplanet(nplanet, jmax, ntrans, params[1:18])
	# for k=1:ntrans[1]
	# 	ttv1[k] = ttv[1,k]
	# end
	for i=1:ntrans[2]
		  #tmax = param[end-2]; phi0 = param[end-1]; deltaphi = param[end]
		  ts = params[end-2] #tmax sinphi0
		  tc = params[end-1] #tmax cosphi0
		  deltaphi = params[end]
		  moon_ttvs = ttv[2,i] + ts*cos((i-1)*deltaphi) + tc*sin((i-1)*deltaphi)
	end
	return moon_ttvs
end