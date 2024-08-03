include("CGS.jl")
using TTVFaster,DataFrames,CSV,LsqFit,Statistics,JLD2,PyPlot

function chisquare(tt0,nplanet,ntrans,params,tt,sigtt,jmax,EM)
  chisq = 0.0  #check memory allocation >>>>>>>>>>>>
  tt_model = ttv_wrapper(tt0,nplanet,ntrans,params,jmax,EM) 
  for j=1:length(tt)
    chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
  end
  return chisq
end
avg(x,y)=(x + y)/2
gaussian(x,mu,sig)=exp.(-((x .- mu).^2) ./ (2 * sig^.2))
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
    @show BIC
    return reduced_chisq, BIC,chisq
  end

function fit_BIC(tt0,nplanet,ntrans,params,tt,sigtt,jmax,EM)
	chi2 = chisquare(tt0,nplanet,ntrans,params,tt,sigtt,jmax,EM)
	N=length(tt0) ; k=length(params)
	BIC=chi2 + k*ln(N)
	return chi2,BIC
end
G=CGS.GRAV /1e3 #in MKS units
AU=CGS.AU /1e2 #in MKS units

Kepler_law(Per,mp,mstar)=((G*(mstar + mp)* (Per*24*3600)^2) /(4*pi^2))^(1/3) 
Hill_radius(Per,mp,ecc,mstar)=(Kepler_law(Per,mp,mstar) * (1-ecc) * (mp/(3 * mstar))^(1/3)) / AU
RV_semiamplitude(Per,mp,ecc,inc,mstar) = (((mp*sin(inc))^3/((Per*24*3600) * (1-ecc^2)^3/2)) * (2pi*G/(mstar+mp)^2))^(1/3)
calc_deg(value)=value * 180/pi
calc_evec1(e,omega)=e* cos(omega-77)
calc_evec2(e,omega)=e* sin(omega-77)
calc_ecc(ecosomega,esinomega)=sqrt(ecosomega^2 + esinomega^2)
calc_tmax(a_p,a_s,m_p,m_s,P_p)=(a_s*m_s*P_p) / (2*pi*a_p*(m_s+m_p))
#Hill_radius((1733*24*3600),(27*5.9742e24),0.4,1.99e30)
function mutual_Hill(Per1,mp1,mstar,Per2,mp2)
	a1,a2 = Kepler_law(Per1,mp1,mstar),Kepler_law(Per2,mp2,mstar)	
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

	# function peaks(df)
  #   for i=1:length(df.lprob)
  #    	if xprob(df.lprob)[i] > .8 && xprob(df.lprob)[i]!=1
  #      println("indx: ",i," Lprob: ", xprob(df.lprob)[i])
  #      append!(indxs,i)
  #     end
  #    end
  # end
  peaks(df)
  new_params=df[indxs,:]
	return new_params
end
function calc_quad_errs(xcos,xcos_err,xsin,xsin_err)
	x = sqrt(xcos^2 .+ xsin^2)
	return sqrt(((xcos^2 * xcos_err^2) + (xsin^2 * xsin_err^2))/x^2)
end

parname=[
	L"$m_b / M_{\odot}$",L"$P_b$",L"$t_{0,b}$",L"$e_b cos(ω_b)$",L"$e_b sin(ω_b)$",
    	L"$m_c / M_{\odot}$",L"$P_c$",L"$t_{0,c}$",L"$e_c cos(ω_c)$",L"$e_c sin(ω_c)$",
 #  	L"$m_e / M_{\odot}$",L"$P_e$",L"$t_{0,e}$",L"$e_e cos(ω_e)$",L"$e_e sin(ω_e)$",
   	L"$m_d / M_{\odot}$",L"$P_d$",L"$t_{0,d}$",L"$e_d cos(ω_d)$",L"$e_d sin(ω_d)$",
#   	L"$μ_5$ ",L"$P_5$ [days]",L"$t_{0,5}$",L"$e_4 cos(ω_5)$",L"$e_5 sin(ω_5)$",
    	L"$t_{max} sin(ϕ_0)$",L"$t_{max} cos(ϕ_0)$",L"$Δϕ$ [rad]",L"$σ_{sys}^2$ [days]"]
truem1,truem2,truem3,truem4=0.815.*CGS.MEARTH/CGS.MSUN,1.0.*CGS.MEARTH/CGS.MSUN,0.1074.*CGS.MEARTH/CGS.MSUN,317.8.*CGS.MEARTH/CGS.MSUN
truep1,truep2,truep3,truep4=224.7007992,365.2564,686.9795859,4332.82012875
trueec1,trueec2,trueec3,trueec4=calc_evec1(0.00677323,131.53298),calc_evec1(0.01671022,102.94719),calc_evec1(0.09341233,336.04084),calc_evec1(0.04839266,14.75385)
truees1,truees2,truees3,truees4=calc_evec2(0.00677323,131.53298),calc_evec2(0.01671022,102.94719),calc_evec2(0.09341233,336.04084),calc_evec2(0.04839266,14.75385)
truee1,truee2,truee3,truee4=0.00677323,0.01671022,0.09341233,0.04839266
true_vals=[truem1;truep1;0.0;trueec1;truees1;truem2;truep2;0.0;trueec2;truees2;
#    truem3;truep3;0.0;trueec3;truees3;
truem4;truep4;0.0;trueec4;truees4]
