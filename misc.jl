include("CGS.jl")
using TTVFaster,DataFrames,CSV,LsqFit
calc_deg(value)=value * 180/pi
calc_rad(value)=value * pi/180
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

