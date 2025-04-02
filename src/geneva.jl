J2000 = 2451545.0
jd1=2.4332825e6
tref=2.43e6
Pers =[224.70 365.25 686.980 4332.82 10755.5]
t0=[3503.7655 3624.4022]
# ecosw=[0.0008 0.013 -0.059 0.034]
# esinw=[-0.002 0.0012 -0.131 -0.033]
omegas=[131.767 102.93 -23.917 14.274 92.861]
eccs=[0.0068 0.0167 0.0934 0.0485 0.0555]
ecoss= eccs.*cosd.(omegas); esins=eccs.*sind.(omegas)
lambdas=[181.979 100.467 -4.553 34.396]
rprstars=[0.00870 0.00916 0.00487 0.10052]


p2= jldopen("MCMC/newp2_mcmc30s28yrs.jld2","r")

iburn=p2["iburn"];nsteps=p2["nsteps"];nwalkers=p2["nwalkers"];nparam=length(p2["pname"])

indepsamples=p2["indepsamples"]

par_mcmc=p2["par_mcmc"][:,iburn:end,:];#(75, 8494, 11)


# for istep = 1:size(par_mcmc)[2]
# trimed_pars=zeros(1950,11)
trimmed_pars = vec([par_mcmc[:,begin:330:end,i]) for i in 1:11]
calc_omega(pomega,Omega) =pomega - Omega
# calc_anom(t,t0,P)=(360 ./ P) .* (t .- t0)
calc_M(t,t0,P)=2pi .* (t.-t0) ./ P # mean anomaly

calc_λ(t,t0,P,esinw)=((360/P) .* (t .- (t0 .+tref)) .+ 2 * esinw # mean longitude
calc_lam(P,t0,t,lambda)=((t0-t)*(360/P)) + lambda
new_lambda(t,P,lambda)=lambda + ((2pi *t)/P)

df=DataFrame(walk_step=collect(1:1950),
		λ_1=mod.(calc_λ.(J2000,trimmed_pars[3],trimmed_pars[2],trimmed_pars[5]),360),
		P_1=trimmed_pars[2],ecos_1=trimmed_pars[4],esin_1=trimmed_pars[5],
		inc1=ones(1950).*90, Ω_1=zeros(1950),mu_1=trimmed_pars[1],
		rprstar_1=ones(1950).*rprstars[1],
	
		λ_2=mod.(calc_λ.(J2000,trimmed_pars[8],trimmed_pars[7],trimmed_pars[10]),360),
		P_2=trimmed_pars[7],ecos_2=trimmed_pars[9],esin_2=trimmed_pars[10],
		inc2=ones(1950).*90,Ω_2=zeros(1950),mu_2=trimmed_pars[6],
		rprstar_2=ones(1950).*rprstars[2],
		mstar=ones(1950),rstar=ones(1950),
# extra column for timing uncertainty
		σ_TT = sqrt.((trimmed_pars[11].*24*3600).^2 .+ 30^2))

samples=string("Lindor_Sol_0_samples.csv")
CSV.write(samples,df)
								# P_3=param_p3[12,:],t03=param_p3[13,:],ecos3=param_p3[14,:],esin3=param_p3[15,:],
								# lprob=lprob_p3[:])