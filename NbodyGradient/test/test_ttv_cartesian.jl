

#include("../src/ttv.jl")
#include("/Users/ericagol/Computer/Julia/regress.jl")

@testset "ttv_cartesian" begin

# This routine takes derivative of transit times with respect
# to the initial cartesian coordinates of bodies. [x]
#n = 8
n = 3
H = [3,1,1]
t0 = 7257.93115525-7300.0
#h  = 0.12
h  = 0.04
#h  = 0.02
#tmax = 600.0
#tmax = 1000.0
#tmax = 100.0
tmax = 10.0

# Read in initial conditions:
elements = readdlm("elements.txt",',')
elements[:,3] .-= 7300.0

# Make an array, tt,  to hold transit times:
# First, though, make sure it is large enough:
ntt = zeros(Int64,n)
for i=2:n
  ntt[i] = ceil(Int64,tmax/elements[i,2])+3
end
println("ntt: ",ntt)
tt  = zeros(n,maximum(ntt))
tt1 = zeros(n,maximum(ntt))
tt2 = zeros(n,maximum(ntt))
tt3 = zeros(n,maximum(ntt))
# Save a counter for the actual number of transit times of each planet:
count = zeros(Int64,n)
count1 = zeros(Int64,n)
# Call the ttv function:
rstar = 1e12
dq = ttv_elements!(H,t0,h,tmax,elements,tt1,count1,0.0,0,0,rstar)
# Now call with half the timestep:
count2 = zeros(Int64,n)
count3 = zeros(Int64,n)
dq = ttv_elements!(H,t0,h/2,tmax,elements,tt2,count2,0.0,0,0,rstar)

# Now, compute derivatives (with respect to initial cartesian positions/masses):
dtdq0 = zeros(n,maximum(ntt),7,n)
dtdelements = ttv_elements!(H,t0,h,tmax,elements,tt,count,dtdq0,rstar)
#read(STDIN,Char)

# Check that this is working properly:
for i=1:n
  for j=1:count2[i]
#    println(i," ",j," ",tt[i,j]," ",tt2[i,j]," ",tt[i,j]-tt2[i,j]," ",tt1[i,j]-tt2[i,j])
  end
end
#read(STDIN,Char)

# Compute derivatives numerically:
# Compute the numerical derivative:
dtdq0_num = zeros(BigFloat,n,maximum(ntt),7,n)
itdq0 = zeros(Int64,n,maximum(ntt),7,n)
dlnq = big(1e-15)
hbig = big(h); t0big = big(t0); tmaxbig=big(tmax); tt2big = big.(tt2); tt3big = big.(tt3)
for jq=1:n
  for iq=1:7
    elements2  = big.(elements)
    dq_plus = ttv_elements!(H,t0big,hbig,tmaxbig,elements2,tt2big,count2,dlnq,iq,jq,big(rstar))
    elements3  = big.(elements)
    dq_minus = ttv_elements!(H,t0big,hbig,tmaxbig,elements3,tt3big,count3,-dlnq,iq,jq,big(rstar))
    for i=1:n
      for k=1:count2[i]
        # Compute double-sided derivative for more accuracy:
        dtdq0_num[i,k,iq,jq] = (tt2big[i,k]-tt3big[i,k])/(dq_plus-dq_minus)
#        println(i," ",k," ",iq," ",jq," ",tt2big[i,k]," ",tt3big[i,k]," ")
      end
    end
  end
end

nbad = 0
ntot = 0
diff_dtdq0 = zeros(n,maximum(ntt),7,n)
mask = zeros(Bool, size(dtdq0))
for i=2:n, j=1:count[i], k=1:7, l=1:n
  if abs(dtdq0[i,j,k,l]-dtdq0_num[i,j,k,l]) > 0.1*abs(dtdq0[i,j,k,l]) && ~(abs(dtdq0[i,j,k,l]) == 0.0  && abs(dtdq0_num[i,j,k,l]) < 1e-3)
#    println(i," ",j," ",k," ",l," ",dtdq0[i,j,k,l]," ",dtdq0_num[i,j,k,l]," ",itdq0[i,j,k,l])
    nbad +=1
  end
  diff_dtdq0[i,j,k,l] = abs(dtdq0[i,j,k,l]-dtdq0_num[i,j,k,l])
  if k != 2 && k != 5
    mask[i,j,k,l] = true
  end
  ntot +=1
end

tt_big = big.(tt); elementsbig = big.(elements); rstarbig = big(rstar)
dqbig = ttv_elements!(H,t0big,hbig,tmaxbig,elementsbig,tt_big,count,big(0.0),0,0,rstarbig)
# Now halve the time steps:
tt_big_half = copy(tt_big)
dqbig = ttv_elements!(H,t0big,hbig/2,tmaxbig,elementsbig,tt_big_half,count1,big(0.0),0,0,rstarbig)

# Compute the derivatives in BigFloat precision to see whether finite difference
# derivatives or Float64 derivatives are imprecise at the start:
dtdq0_big = zeros(BigFloat,n,maximum(ntt),7,n)
hbig = big(h); tt_big = big.(tt); elementsbig = big.(elements); rstarbig = big(rstar)
dtdelements_big = ttv_elements!(H,t0big,hbig,tmaxbig,elementsbig,tt_big,count,dtdq0_big,rstarbig)

#=
using PyPlot

clf()
# Plot the difference in the TTVs:
for i=2:3
#  diff1 = abs.(tt1[i,2:count1[i]]./tt_big[i,2:count1[i]]-1.0);
  diff1 = convert(Array{Float64,1},abs.(tt1[i,2:count1[i]].-tt_big[i,2:count1[i]])/elements[i,2]);
  dtt=tt[i,2:count1[i]].-tt[i,1]
  loglog(dtt,diff1);
#  diff2 = abs.(tt2[i,2:count1[i]]./tt_big_half[i,2:count1[i]]-1.0);
#  diff2 = abs.(tt2[i,2:count1[i]].-tt_big_half[i,2:count1[i]])/elements[i,2];
#  loglog(tt[i,2:count[i]]-tt[i,1],diff2);
end
loglog([1.0,1024.0],2e-15*[1,2^15],":")
for i=2:3, k=1:7, l=1:3
  if maximum(abs.(dtdelements_big[i,2:count[i],k,l])) > 0
    diff1 = convert(Array{Float64,1},abs.(dtdelements_big[i,2:count[i],k,l]./dtdelements[i,2:count[i],k,l].-1));
    diff3 = convert(Array{Float64,1},abs.(asinh.(dtdelements_big[i,2:count[i],k,l])-asinh.(dtdelements[i,2:count[i],k,l])));
#    loglog(tt[i,2:count[i]]-tt[i,1],diff1);
    dtt = tt[i,2:count[i]].-tt[i,1]
    loglog(dtt,diff3);
    println(i," ",k," ",l," frac error: ",convert(Float64,maximum(diff1))," asinh error: ",convert(Float64,maximum(diff3))); #read(STDIN,Char);
  end
  if maximum(abs.(dtdq0_big[i,2:count[i],k,l])) > 0
    diff1 = convert(Array{Float64,1},abs.(dtdq0[i,2:count[i],k,l]./dtdq0_big[i,2:count[i],k,l].-1.0));
#    diff2 = abs.(asinh.(dtdq0_big[i,2:count[i],k,l]).-asinh.(dtdq0_num[i,2:count[i],k,l]));
    diff3 = convert(Array{Float64,1},abs.(asinh.(dtdq0_big[i,2:count[i],k,l]).-asinh.(dtdq0[i,2:count[i],k,l])));
#    loglog(tt[i,2:count[i]].-tt[i,1],diff1);
    dtt = tt[i,2:count[i]].-tt[i,1]
    loglog(dtt,diff3,linestyle=":");
#    loglog(tt[i,2:count[i]].-tt[i,1],diff2,".");
    println(i," ",k," ",l," frac error: ",convert(Float64,maximum(diff1))," asinh error: ",convert(Float64,maximum(diff3))); #read(STDIN,Char);
  end
end
mederror = zeros(size(tt))
for i=2:3
  for j=1:count1[i]
    data_list = Float64[]
    for k=1:7, l=1:3
      if abs(dtdq0_num[i,j,k,l]) > 0
        push!(data_list,abs(dtdq0[i,j,k,l]/dtdq0_num[i,j,k,l]-1.0))
      end
    end
    mederror[i,j] = median(data_list)
  end
end

# Plot a line that scales as time^{3/2}:

loglog([1.0,1024.0],1e-12*[1,2^15],":",linewidth=3)

=#
#println("Max diff log(dtdq0): ",maximum(abs.(dtdq0_num[mask]./dtdq0[mask].-1.0)))
println("Max diff asinh(dtdq0): ",maximum(abs.(asinh.(dtdq0_num[mask]).-asinh.(dtdq0[mask]))))
#unit = ones(dtdq0[mask])
#@test isapprox(dtdq0[mask]./convert(Array{Float64,4},dtdq0_num)[mask],unit;norm=maxabs)
#@test isapprox(dtdq0[mask],convert(Array{Float64,4},dtdq0_num)[mask];norm=maxabs)
@test isapprox(asinh.(dtdq0[mask]),asinh.(convert(Array{Float64,4},dtdq0_num)[mask]);norm=maxabs)
#end

#
#nderiv = n^2*7*maximum(ntt)
#loglog(abs.(reshape(dtdq0,nderiv)),abs.(reshape(convert(Array{Float64,4},dtdq0_num),nderiv)),".")
#loglog(abs.(reshape(dtdq0,nderiv)),abs.(reshape(convert(Array{Float64,4},diff_dtdq0),nderiv)),".")


## Make a plot of some TTVs:
#
#fig,axes = subplots(4,2)
#
#for i=2:8
#  ax = axes[i-1]
#  fn = zeros(Float64,2,count1[i])
#  sig = ones(count1[i])
#  tti1 = tt1[i,1:count1[i]]
#  tti2 = tt2[i,1:count2[i]]
#  for j=1:count1[i]
#    fn[1,j] = 1.0
#    fn[2,j] = round(Int64,(tti1[j]-elements[i,3])/elements[i,2])
#  end
#  coeff,cov = regress(fn,tti1,sig)
#  tt_ref1 = coeff[1]+coeff[2]*fn[2,:]
#  ttv1 = (tti1-tt_ref1)*24.*60.
#  coeff,cov = regress(fn,tti2,sig)
#  tt_ref2 = coeff[1]+coeff[2]*fn[2,:]
#  ttv2 = (tti2-tt_ref2)*24.*60.
#  ax[:plot](tti1,ttv1)
##  ax[:plot](tti2,ttv2)
#  ax[:plot](tti2,((ttv1-ttv2)-mean(ttv1-ttv2)))
#  println(i," ",coeff," ",elements[i,2:3]," ",coeff[1]-elements[i,3]," ",coeff[2]-elements[i,2])
#  println(i," ",maximum(ttv1-ttv2-mean(ttv1-ttv2))*60.," sec ", minimum(ttv1-ttv2-mean(ttv1-ttv2))*60.," sec" )
end
