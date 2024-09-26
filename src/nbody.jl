using Photodynamics,PyPlot,Random

module AgolModels # will combine with common.jl
   using Photodynamics
   # struct AbstractModel <: AbstractFloat end
   mutable struct PhotometryModel{T<:Real} #<: AbstractModel
      θ::Vector{T}
      ic::ElementsIC{T}
      tt::TransitTiming{T}
      ts::TransitSeries{T,Photodynamics.ComputedTimes}
      d::NbodyGradient.Derivatives{T}
      lc::Lightcurve{T}
      intr::Integrator
      npar_per_pl::Int64
      N::Int64
      t0::T
      tmax::T
      h::T
   end
   # do_grad::Bool
   #     J::Matrix{T}
   """ 
   We can also make immutable PhotometryModel structure that doens't have the θ vector.
      θ=[pl_params;k;u_n;rstar;mstar]
      - k  vector has radius ratios
      - u_n vector has limb dark coefficients
   """
   function PhotometryModel(θ,nplanet,t0,tmax,h;vary_orb_angles::Bool=false) 
       if vary_orb_angles; npar_per_pl=7 ; else npar_per_pl=5 end
       N=nplanet+1
       @assert(length(θ)> npar_per_pl*nplanet)
       k = sqrt.(θ[npar_per_pl*nplanet+1:end-4]) #radius ratios
       u_n = (θ[end-3:end-2]) # quadratic limb dark. coeefs
       rstar = θ[end-1]
       mstar = θ[end] 
       elements=zeros(nplanet+1,7)
       elements[1,1]=1.0 # assume that were using mass-ratios, but just in case
       for i=1:nplanet 
           elements[i+1,1:5] .= θ[(i-1)*5 + 1:5*i]
           if vary_orb_angles 
            @assert(npar_per_pl==7) #  Inclination and Ω must be provided
           elements[i+1,6:7] =  θ[(i-1)*5+6:(i-1)*5+7]
           else
            # Hold the inclination fixed at 90° and the Ω at π
           elements[i+1,6] = pi/2;    
           elements[i+1,7] = pi 
           end
       end
       lc = Lightcurve(h, tmax-t0, copy(u_n), copy(k), rstar);
       ic = ElementsIC(t0,N,elements)
       intr=Integrator(h,t0,tmax)
       tt = TransitTiming(intr.tmax, ic)
       ts = TransitSeries(intr.tmax, ic)
       d = NbodyGradient.Derivatives(Float64, ic.nbody);
    return PhotometryModel(θ,ic,tt,ts,d,lc,intr,npar_per_pl,N,t0,tmax,h)
   end

   export PhotometryModel
end
H=[-1 1 0 0 0;0 0 -1 1 0; -1 -1 1 1 0; -1 -1 -1 -1 1;-1 -1 -1 -1 -1] #heirarchy matrix

star=Elements(m=1.0)
p1=Elements(m=2.5422162292092204e-6,
          P=224.70078014619864,
         t0=3503.765349062694,
           ecosω=-0.003,
           esinω=-0.006,
           I=π/2,
           Ω=0.0);
p2=Elements(m=3.0256056455411807e-6,
          P=365.2564540136157,
         t0=3624.4021734758985,
            ecosω=0.011,
            esinω=0.012,
               I=π/2,
               Ω=0.0);
# p4=Elements(m=0.3227e-6,P= 686.980 ,t0=383.823 , ecosω =0.0403  ,  esinω= 0.0268 );

moon=Elements(m=0.012*(p2.m),P=27.3,t0=3624.3949611424468,
e=0.055,I=(pi/2)-deg2rad(5),Ω=deg2rad(60))
# wont get good constraint on moon's orbital elements.
# actual inc is 5deg wrt ecliptic. could be pos or neg
# S-V ; E-M ; SV_bc-EM_bc; must define inner binaries before outer binaries


t0=3283.5;tmax=4000.0;
ic=ElementsIC(t0,[star;p1;p2;moon;p4])
rstar = 0.00465047 # Sol (Rstar in AU)
u_n = [0.39256, 0.29064];  # Quad. Limbdarkening coefficients
k_moon = [0.0087002, 0.0091705, 0.27*0.0091705]  # Radius ratios with moon

cadence = 2 / 60 / 24  # 2 minute cadence in days
obs_duration = tmax-t0  # Duration of observations in days
tobs = collect(t0:cadence:t0+obs_duration) 
# allocate arrays and initialzie state
lc = Lightcurve(cadence, tobs, ones(length(tobs)), zeros(length(tobs)), u_n, k_moon, rstar); 
intr = Integrator(0.05, t0, tmax) ;
s=State(ic)
ts = TransitSeries(intr.tmax,ic) ;# args: max time obs, IC, h [step size for series pnts]
tt = TransitTiming(intr.tmax,ic) ;#TransitTiming(intr.tmax,ic) # ;need this for grad;
intr(s, ts, tt; grad=false);
compute_lightcurve!(lc, ts)
# plot(lc.tobs,lc.flux)
d = NbodyGradient.Derivatives(Float64, ic.nbody);
# @show tt
# model=PhotometryModel(ic,tt,ts,d,lc);
Random.seed!(42)
dy=1 .+ randn(length(tobs))
fig=figure(figsize=(7,5));ax=fig.add_subplot(111);
ax.plot(lc.tobs,lc.flux);
# ax.errorbar(lc.tobs,lc.flux,yerr=dy,fmt="."); # hang-up 
ax.set_ylabel("Relative Flux");ax.set_xlabel("Time [days]")
fig.savefig("testing.png")
