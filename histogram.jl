
function histogram(param,nbin)
  p1 = minimum(param)-1e-15; p2 = maximum(param)+1e-15
  pbin_square = [p1]
  hist_square = [0.0]
  pbin = zeros(nbin)
  hist = zeros(nbin)
  psort = sort(param,dims=1)
  i1 = 1; np = size(param)[1]
  #println(p1," ",psort[i1]," ",p2," ",psort[np])
  for i=1:nbin
    while psort[i1] <= p1+(p2-p1)/nbin*i
      hist[i] += 1.0
      i1 += 1
      if i1 == np
        break
      end
    end
    push!(hist_square,hist[i]); push!(hist_square,hist[i])
    pbin[i] = p1+(p2-p1)/nbin*(i-0.5)
    push!(pbin_square,p1+(p2-p1)/nbin*(i-1))
    push!(pbin_square,p1+(p2-p1)/nbin*i)
    if i1 == np
      break
    end
  end
  push!(hist_square,0.0)
  push!(pbin_square,p2)
  return pbin,hist,pbin_square,hist_square
end

# function plot_hist(param,nbin)

# end
# x = randn(100000)
# xbin,xhist,xbin_square,xhist_square = histogram(x,50)
# plot(xbin,xhist)
# plot(xbin_square,xhist_square)