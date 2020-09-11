# Julia v1.1
include("ttv_nplanet.jl")
function ttv_wrapper(tt0, nplanet, ntrans, params, fixp3::Bool = false, p3_cur::Float64 = 0.0, EMB::Bool = true)
    # These lines need modification for different choices of parameters:
    if nplanet == 2
        n1, n2 = ntrans
    end
    if nplanet == 3
        n1, n2, n3 = ntrans
    end
    if nplanet == 3 && fixp3
        param = [params[1:11];p3_cur;params[12:end]]
    else 
        param = params
    end
    jmax = 5
    # Call ttv_nplanet:
    ttv = ttv_nplanet(nplanet, jmax, ntrans, param[1:5*nplanet])
    # We measure transit times, not TTVs, so add back in the linear ephemeris:
    t01 = params[3]
    per1 = params[2]
    ttv1 = collect(range(t01,stop = t01+per1*(n1-1),length = n1)) #this doesnt account for skipped transits
    for i=1:n1
     ttv1[i]+= ttv[1,i]
    end
    t02 = params[8]
    per2 = params[7]
    ttv2 = collect(range(t02,stop = t02+per2*(n2-1),length = n2))
    for i=1:n2
        if EMB
            ttv2[i] += ttv[2,i]
        else
            #tmax = param[end-2]; phi0 = param[end-1]; deltaphi = param[end]
            ts = param[end-2] #tmax sinphi0
            tc = param[end-1] #tmax cosphi0
            deltaphi = param[end]
            ttv2[i] += ttv[2,i] + ts*cos((i-1)*deltaphi) + tc*sin((i-1)*deltaphi)
        end
    end
    # If transit times of additional planets were observable these would need to be added in.
    #println("param2: ",param)
    return [ttv1;ttv2]
end


