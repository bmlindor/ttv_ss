# Computes TTVs with TTVFaster for N planets with pairwise TTV calculation.
include("compute_ttv.jl")

function ttv_nplanet(nplanet::Int64,jmax::Int64,ntrans::Vector{Int64},data::Vector,sqrte::Bool)
  # Need at least two planets!
  @assert(nplanet>=2)
  # The ntrans vectors should have length nplanet:
  @assert(length(ntrans)==nplanet)
  # Define type of ttv array:
  ttv_el_type = eltype(data) == Float64 ? Float64 : Number
  # Need to create an array to store TTVs with maximum length equal to maximum number
  # of transit times of any planet:
  ntransmax = maximum(ntrans)
  #ttv = zeros(ttv_el_type,nplanet,ntransmax)
  ttv = zeros(nplanet,ntransmax)
  # Each planet requires 5 elements in data: mass_ratio,period,trans0,ecosw,esinw:
  @assert(length(data)==5*nplanet)
  @assert(jmax>=1)  # Should there be a larger minimum?
  for iplanet=1:nplanet
  # Each planet should have at least 2 transits:
    @assert(ntrans[iplanet]>=2)
  end
  for iplanet=1:nplanet-1
  # The periods of the planets should be ordered from least to greatest:
    @assert(data[(iplanet-1)*5+2] < data[iplanet*5+2])
  end
  # Set up planets planar-planet types for all of the planets:
  #planet = Array{Planet_plane_hk}(nplanet)
  #planet = Array{Any}(nplanet)
  # Loop over pairs of planets to compute pairwise TTVs
  # Loop over inner planets:
  #println("Looping over planets in ttv_nplanet:")
  for iplanet=1:nplanet-1
    # Create a Planet_plane_hk type for the inner planet:
    if sqrte
      p1=TTVFaster.Planet_plane(data[(iplanet-1)*5+1],data[(iplanet-1)*5+2],data[(iplanet-1)*5+3],data[(iplanet-1)*5+4],data[(iplanet-1)*5+5])
    else
      p1=TTVFaster.Planet_plane_hk(data[(iplanet-1)*5+1],data[(iplanet-1)*5+2],data[(iplanet-1)*5+3],data[(iplanet-1)*5+4],data[(iplanet-1)*5+5])
    end
    # Create an array of times for the inner planet:
    n1 = ntrans[iplanet]
    time1 = collect(p1.trans0 .+ range(0,stop=n1-1,length=n1) .* p1.period)
    # Loop over outer planets:
    for jplanet=iplanet+1:nplanet
      # Create a Planet_plane_hk type for the outer planet:
      if sqrte
        p2=TTVFaster.Planet_plane(data[(jplanet-1)*5+1],data[(jplanet-1)*5+2],data[(jplanet-1)*5+3],data[(jplanet-1)*5+4],data[(jplanet-1)*5+5])
      else
        p2=TTVFaster.Planet_plane_hk(data[(jplanet-1)*5+1],data[(jplanet-1)*5+2],data[(jplanet-1)*5+3],data[(jplanet-1)*5+4],data[(jplanet-1)*5+5])
      end 
      # Create an array of times for the outer planet:
      n2 = ntrans[jplanet]
      time2 = collect(p2.trans0 .+ range(0,stop=n2-1,length=n2) .* p2.period)
      # Define arrays to hold the TTVs:
      ttv1=zeros(ttv_el_type,n1)
      ttv2=zeros(ttv_el_type,n2)
      # Call the compute_ttv code which implements equation (33) from Agol & Deck (2016):
      #    println("Calling compute_ttv")
      TTVFaster.compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2)
      #    println("Finished compute_ttv")
      for i=1:n1
        ttv[iplanet,i] += ttv1[i]
      end
      for i=1:n2
        ttv[jplanet,i] += ttv2[i]
      end
    end
  end
  return ttv
end

function ttv_wrapper(tt0,nplanet,ntrans,params,jmax,sqrte::Bool=false,EMB::Bool=true)
  # These lines need modification for different choices of parameters:
  if nplanet == 2
    n1,n2 = ntrans
  end
  if nplanet == 3
    n1,n2,n3 = ntrans
  end
  # if nplanet == 3 && fixp3
  #     param = [params[1:11];p3_cur;params[12:end]]
  # else 
  #     param = params
  # end
  # jmax = 5
  # Call ttv_nplanet:
  ttv = ttv_nplanet(nplanet,jmax,ntrans,params[1:5*nplanet],sqrte)
  # We measure transit times,not TTVs,so add back in the linear ephemeris:
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
      ts = params[end-2] #tmax sinphi0
      tc = params[end-1] #tmax cosphi0
      deltaphi = params[end]
      ttv2[i] += ttv[2,i] + ts*cos((i-1)*deltaphi) + tc*sin((i-1)*deltaphi)
    end
  end
  # If transit times of additional planets were observable these would need to be added in.
  #println("param2: ",param)
  return [ttv1;ttv2]
end

function chisquare(tt0,nplanet,ntrans,params,tt,sigtt,jmax,sqrte)#,fixp3::Bool = false,p3_cur::Float64 = 0.0)
  chisq = 0.0
  # println(params,tt[1],sigtt[1])
  tt_model = ttv_wrapper(tt0,nplanet,ntrans,params,jmax,sqrte) #,fixp3,p3_cur)
  for j=1:length(tt)
    chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
  end
  # println(nplanet)
  return chisq
end