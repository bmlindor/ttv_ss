using Photodynamics,PyPlot,Random

struct PhotometryModel{T<:Real}
           ic::ElementsIC{T}
           tt::TransitTiming{T}
           ts::TransitSeries{T,Photodynamics.ComputedTimes}
           d::NbodyGradient.Derivatives{T}
           lc::Lightcurve{T}
       #     intr::Integrator{T}
       #     J::Matrix{T}
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

moon=Elements(m=0.012*(p2.m),P=27.3,t0=3624.3949611424468,ecosω=p2.ecosω,esinω=p2.ecosω,I=π/2);


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
# # model=PhotometryModel(ic,tt,ts,d,lc);
Random.seed!(42)
dy=1 .+ randn(length(tobs))
fig=figure(figsize=(7,5));ax=fig.add_subplot(111);
ax.plot(lc.tobs,lc.flux);
# ax.errorbar(lc.tobs,lc.flux,yerr=dy,fmt=".");
ax.set_ylabel("Relative Flux");ax.set_xlabel("Time [days]")
fig.savefig("../IMAGES/testing.png")
