
# This code tests the autodiffed Drift/Kepler step

include("../src/kepler_drift_gamma.jl")

@testset "kepler_drift_gamma" begin


function test_kepler_drift(dlnq::BigFloat,drift_first::Bool)
# Call as: save,jac_num,jacobian=test_kepler_drift(1e-4,true)
# This routine runs a test of the kep_elliptic_jacobian function in kepler_solver.jl
# Define the central force constant in terms of AU and days:
k = (2.0*pi/365.25)^2
# Initial position at 1 AU:
x0 = [0.01,1.0,0.01]
# Circular velocity:
r0 = sqrt(x0[1]*x0[1]+x0[2]*x0[2]+x0[3]*x0[3])
vcirc = sqrt(k/r0)
# Define initial velocity at apastron:
v0 = [.9*vcirc,0.01*vcirc,0.01*vcirc]  # The eccentricity is about ~2(1-v0/vcirc).
#h = 18.0 # 18-day timesteps
h = 0.000005 # 18-day timesteps
#h = 0.05 # 18-day timesteps
hbig = big(h)

s0::Float64 = 0.0
x = zeros(Float64,3)
v = zeros(Float64,3)
t = 0.0
nsteps = 1
xsave=zeros(Float64,12,nsteps)  # First is time (1,:); next three are position (2-4,:); next three velocity (5-7,:); then r (8,:), drdt (9,:); then beta (10,:); finally s & ds (11:12,:)
# Create a variable to store the state at the end of a step:
state=zeros(Float64,12)
state_diffbig=zeros(BigFloat,12)

xsave[1,1]=0.0
for jj=1:3
  xsave[jj+1,1]=x0[jj]
  xsave[jj+4,1]=v0[jj]
end
# Save beta:
xsave[8,1]=r0
beta0 = 2.0*k/r0-dot(v0,v0)
xsave[10,1]=beta0
# First, compute it with the old version:
jacobian_old =zeros(Float64,7,7)
d_old = Derivatives(Float64)
iter = kep_drift_ell_hyp!(x0,v0,k,h,s0,state,jacobian_old,drift_first,d_old)
#println("s old: ",state[11])

# Check that old and new are giving the same answer:
delxv = jac_delxv_gamma!(x0,v0,k,h,drift_first;debug=true)
println("state: ",state)
println("delxv: ",delxv)
println("Old-new: ",state[2:7]-delxv[1:6])

# Now compute algebraic Jacobian:
delxv,jac_alg = jac_delxv_gamma!(x0,v0,k,h,drift_first;grad=true,auto=false,debug=true)

# Next, compute autodiff Jacobian:
delxv,jacobian = jac_delxv_gamma!(x0,v0,k,h,drift_first;grad=true,auto=true,debug=true)


# Now, do finite differences at higher precision:
kbig = big(k)
x0big = big.(x0); v0big = big.(v0)
# Now compute big-float precision autodiff Jacobian:
delxv_big,jacobian_big = jac_delxv_gamma!(x0big,v0big,kbig,hbig,drift_first;grad=true,auto=true,debug=true)
#println("Auto diff Jacobian: ",jacobian)
jac_frac = jac_alg[1:6,:]./convert(Array{Float64,2},jacobian_big[1:6,:]).-1.0
println("Fractional Jacobian difference: ",maxabs(jac_frac[.~isnan.(jac_frac)]))


# Now compute finite-difference Jacobian:
delxv,jac_num = jac_delxv_gamma!(x0big,v0big,kbig,hbig,drift_first;grad=true,auto=false,dlnq=dlnq,debug=true)
name = ["x01","x02","x03","v01","v02","v03","gamma","r","fm1","fdot","gmhgdot","gdotm1"]
for i=1:12
  println(i," ",name[i])
  println("Fini diff Jacobian: ",convert(Array{Float64,1},jac_num[i,:]))
  println("Algebraic Jacobian: ",convert(Array{Float64,1},jac_alg[i,:]))
  println("Auto diff Jacobian: ",convert(Array{Float64,1},jacobian_big[i,:]))
  println("Frac diff Jacobian: ",jac_alg[i,:]./convert(Array{Float64,1},jacobian_big[i,:]).-1.0)
end
println("gamma alg  gradient: ",jac_alg[7,:])
println("gamma auto gradient: ",convert(Array{Float64,1},jacobian_big[7,:]))
println("gamma diff gradient: ",convert(Array{Float64,1},jac_num[7,:]))
println("r alg  gradient: ",jac_alg[8,:])
println("r auto gradient: ",convert(Array{Float64,1},jacobian_big[8,:]))
println("r diff gradient: ",convert(Array{Float64,1},jac_num[8,:]))
println("fm1 gradient  alg: ",jac_alg[9,:])
println("fm1 gradient auto: ",convert(Array{Float64,1},jacobian_big[9,:]))
println("fm1 gradient diff: ",convert(Array{Float64,1},jac_num[9,:]))
println("fdot gradient  alg: ",jac_alg[10,:])
println("fdot gradient auto: ",convert(Array{Float64,1},jacobian_big[10,:]))
println("fdot gradient diff: ",convert(Array{Float64,1},jac_num[10,:]))
println("g-h gdot gradient   alg: ",jac_alg[11,:])
println("g-h gdot gradient  auto: ",convert(Array{Float64,1},jacobian_big[11,:]))
println("g-h gdot gradient  diff: ",convert(Array{Float64,1},jac_num[11,:]))
println("gdot-1 gradient  alg: ",jac_alg[12,:])
println("gdot-1 gradient auto: ",convert(Array{Float64,1},jacobian_big[12,:]))
println("gdot-1 gradient diff: ",convert(Array{Float64,1},jac_num[12,:]))

println("Dim of jacobian_big: ",size(jacobian_big)," dim of jac_num: ",size(jac_num))
println("Maximum jac_big-jac_num: ",maxabs(convert(Array{Float64,2},jacobian_big-jac_num)))
println("Maximum jac_big-jac_alg: ",maxabs(convert(Array{Float64,2},jacobian_big-jac_alg)))
return xsave,jac_num,jac_alg
end

# First try:
xsave,jac_num1,jacobian=test_kepler_drift(big(1e-20),true)
emax = 0.0; imax = 0; jmax = 0
#for i=1:6, j=1:8
for i=1:12, j=1:8
  if jacobian[i,j] != 0.0
    diff = abs(convert(Float64,jac_num1[i,j])/jacobian[i,j]-1.0)
    if  diff > emax
      emax = diff; imax = i; jmax = j
#      println("i ",i," j ",j," jac_num1: ",convert(Float64,jac_num1[i,j])," jacobian: ",jacobian[i,j])
    end
  end
end
println("Maximum fractional error: ",emax," ",imax," ",jmax)
println("Maximum jacobian-jac_num: ",maxabs(jacobian-convert(Array{Float64,2},jac_num1)))

# Second try:
xsave,jac_num1,jacobian=test_kepler_drift(big(1e-20),false)

emax = 0.0; imax = 0; jmax = 0
for i=1:6, j=1:8
  if jacobian[i,j] != 0.0
    diff = abs(convert(Float64,jac_num1[i,j])/jacobian[i,j]-1.0)
    if  diff > emax
      emax = diff; imax = i; jmax = j
#      println("i ",i," j ",j," jac_num1: ",convert(Float64,jac_num1[i,j])," jacobian: ",jacobian[i,j])
    end
  end
end
println("Maximum fractional error: ",emax," ",imax," ",jmax)
println("Maximum jacobian-jac_num: ",maxabs(jacobian-convert(Array{Float64,2},jac_num1)))

@test isapprox(jacobian,jac_num1;norm=maxabs)
end
