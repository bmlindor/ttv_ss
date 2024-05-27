include("CGS.jl")
using TTVFaster,DataFrames,CSV,LsqFit
# calc_omega(pomega,Omega) = pomega - Omega
calc_MeanAnom(t,t0,P)=2pi .* (t.-t0) ./ P
calc_Long(t,t0,P,esinw)=((360/P) .* (t.-t0)) .+ 2*esinw ;
calc_Long(t,t0,P,e,w)=((360/P) .* (t.-t0)) .+ 2*e*sind(w) ;
calc_evec1(e,omega)=e* cos(omega-77)
calc_evec2(e,omega)=e* sin(omega-77)
calc_ecc(ecosom,esinom)=sqrt(ecosom^2 + esinom^2)
calc_tmax(a_p,a_s,m_p,m_s,P_p)=(a_s*m_s*P_p) / (2*pi*a_p*(m_s+m_p))
function chisquare(tt0,nplanet,ntrans,params,tt,sigtt,jmax,EM)
  chisq = 0.0  #check memory allocation >>>>>>>>>>>>
  tt_model = ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM) 
  for j=1:length(tt)
    chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
  end
  return chisq
end
# avg(x,y)=(x + y)/2
gaussian(x,mu,sig)=exp.(-((x .- mu).^2) ./ (2 * sig^.2))
function fit_BIC(tt0,nplanet,ntrans,params,tt,sigtt,jmax,EM)
	chi2 = chisquare(tt0,nplanet,ntrans,params,tt,sigtt,jmax,EM)
	N=length(tt0) ; k=length(params)
	BIC=chi2 + k*ln(N)
	return chi2,BIC
end
  function calc_BIC(lprob,tt0,tt,sigtt,nplanet,ntrans,par_mcmc;EM=false)
    imax=argmax(lprob)
    prob_max=exp.(lprob[imax])
    function calc_chisq(par_mcmc,nplanet,ntrans)
    chisq = 0.0  
    jmax=5
    tt_model = TTVFaster.ttv_wrapper(tt0,nplanet,ntrans,vec(par_mcmc[imax,1:(nplanet*5)+1]),jmax,EM) 
      for j=1:length(tt)
        chisq += (tt[j]-tt_model[j])^2 / (sigtt[j]^2 + par_mcmc[imax,(nplanet*5)+1]^2)
      end
    return chisq
    end
    chisq=calc_chisq(par_mcmc,nplanet,ntrans)
    N=length(tt0) ; k=nplanet*5 + 1
    #println("[N_obs]= ",N," [no. of model params]= ",k)
    #println("chi^2=",chisq)
    # println("max Prob=",prob_max)
    reduced_chisq=chisq/(N-k)
    BIC_chi(chisq,k,N)=chisq + k*log(N)
    BIC=-2*log(prob_max) + k*log(N)
    return reduced_chisq, BIC,chisq
  end
G=CGS.GRAV /1e3 #in MKS units
AU=CGS.AU /1e2 #in MKS units

calc_sma(Per,mp,mstar)=((G*(mstar + mp)* (Per*24*3600)^2) /(4*pi^2))^(1/3) 
Hill_radius(Per,mp,ecc,mstar)=(calc_sma(Per,mp,mstar) * (1-ecc) * (mp/(3 * mstar))^(1/3)) / AU
#Hill_radius((1733*24*3600),(27*5.9742e24),0.4,1.99e30)
RV_semiamplitude(Per,mp,ecc,inc,mstar) = (((mp*sin(inc))^3/((Per*24*3600) * (1-ecc^2)^3/2)) * (2pi*G/(mstar+mp)^2))^(1/3)

function calc_sma(Per::Real,mu::Real) 
  # mp=mu .* (CGS.MSUN/CGS.MEARTH)/CGS.KILOGRAM
  # mstar=CGS.MSUN/CGS.KILOGRAM
  a=((((CGS.GRAV*CGS.MSUN) .+ mu).* (Per.*24*3600).^2) /(4*pi^2)).^(1/3) 
  return a #(AU/1e2) # returns semi-major axis in AU
end
calc_P(a,mu)=
function mutual_Hill(Per1,mu1,Per2,mu2) 
#Hill_stability(μ1::Real,μ2::Real,a1::Real,a2::Real,Cx::Real)
	@assert(Per1 <= Per2)
	a1=calc_sma(Per1,mu1)
	a2=calc_sma(Per2,mu2)
	return ((mp1 + mp2)/(3 * mstar))^(1/3) * (a1 + a2)/2
end
  # Check stabilities
  # stable=true
  # radius=Hill_radius(avg[12],(avg[11].*CGS.MSUN/CGS.MEARTH),calc_ecc(avg[14],avg[15]),CGS.MSUN)
  # ab_mutual_radius=mutual_Hill(avg[2],avg[1].*CGS.MSUN/CGS.MEARTH,CGS.MSUN,avg[7],avg[6].*CGS.MSUN/CGS.MEARTH)
  # bc_mutual_radius=mutual_Hill(avg[7],avg[6].*CGS.MSUN/CGS.MEARTH,CGS.MSUN,avg[12],avg[11].*CGS.MSUN/CGS.MEARTH)
  # println("Planet a-b mutual Hill_radius: ",ab_mutual_radius)
  # println("Planet b-c mutual Hill_radius: ",bc_mutual_radius)
  # println("Planet c Hill radius: ",radius)

function second_peak_params(grid_file::String)
	#Read in as delm
	# data, header=readdlm(grid_file,',',header=true)
	# pnts=length(data[:,1]) ; nparam=length(data[1,:])-1
	# lprob_best=maximum(data[:,end])
	#Read in as data frame (takes longer but better for keeping track)
	df=CSV.read(grid_file,header=1,delim=',',types=Float64,DataFrame)
	max_lprob=maximum(df.lprob)
	# new_params=zeros(nparam)
	prob=exp.(df.lprob .- maximum(df.lprob))
	println("lprob max: ",max_lprob)
	indxs=[]
	function peaks(df)
    for i=1:length(df.lprob)
     	if prob[i] > .8 && prob[i]!=1
       println("indx: ",i," Lprob: ", prob[i])
       append!(indxs,i)
      end
     end
  end
  peaks(df)
  new_params=df[indxs,:]
	return new_params
end
function calc_quad_errs(xcos,xcos_err,xsin,xsin_err)
	x = sqrt(xcos^2 .+ xsin^2)
	return sqrt(((xcos^2 * xcos_err^2) + (xsin^2 * xsin_err^2))/x^2)
end
"""
secular perturbations from Solar System Dynamics book, Appendix A.3 & A.4 has approx rates of change for solar system object per century. Check that these agree with the longitude difference we found of 77 degrees
# elements .+ el_rates.*T
"""
# J2000 = 2451545.0
# T=(jd1-J2000)/36525
# e, I, L ,\vapi, \Omega
elements=[
    0.00677672 3.39467605 181.97909950 131.60246718 76.67984255;
0.01671123 -0.00001531 100.46457166 102.93768193 0.0;
0.09339410 1.84969142 -4.55343205 -23.94362959 49.55953891;
0.04838624 1.30439695 34.39644051 14.72847983 100.47390909;
0.05386179 2.48599187 49.95424423 92.59887831 113.66242448]
el_rates=[
    -0.00004107 -0.00078890 58517.81538729 0.00268329 -0.27769418;
    -0.00004392 -0.01294668 35999.37244981 0.32327364 0.0;
    0.00007882 -0.00813131 19140.30268499 0.44441088 -0.29257343;
    -0.00013253 -0.00183714 3034.74612775 0.21252668 0.20469106;
    -0.00050991 0.00193609 1222.49362201 -0.41897216 -0.28867794
]
