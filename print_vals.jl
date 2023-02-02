using PyPlot,JLD2,Statistics,LsqFit,TTVFaster
# import Main.TTVFaster.ttv_wrapper
include("CGS.jl")
include("misc.jl")
# Print results from MCMC run after burn-in
function print_vals(sigma::Real,nyear::Real,grid_type_nplanet::String,nplanet::Int64,case_num::Int64)
    mcfile=string("MCMC/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    if case_num==1 && isfile(string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
        mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
        case_label="Case 1"
    elseif case_num==2 && isfile(string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
        mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
        case_label="Case 2"
    else
        return println("MCMC file for case ",case_num," with ",grid_type_nplanet," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
    end
    jldmc=jldopen(String(mcfile),"r")
    par_mcmc,lprob_mcmc=f["par_mcmc"],f["lprob_mcmc"]
    iburn,samples=jldmc["iburn"], jldmc["indepsamples"]
    nwalkers,nsteps=jldmc["nwalkers"],jldmc["nsteps"]
    param=f["param"]
    pname=Array{String,length(param)}
    # for i=1:length(param)
    #     if mod(i,5)
    #     pname[i] =
    # prefix=["mu_","P_","t0","ecos","esin",
    #       "mu_2","P_2","t02","ecos2","esin2",
    #       "mu_3","P_3","t03","ecos3","esin3", 
    #       "mu_4","P_4","t04","ecos4","esin4", 
    #       "mu_5","P_5","t05","ecos5","esin5"]
    if model=="p3"
        pname=pname[1:15]
    elseif model=="p4"
        pname=pname[1:20]
    elseif model=="p5"
        pname=pname[1:25]
    end
    # if n
    #   moon_names=["tcosϕ","tsinϕ","Δϕ"]
    #   append!(pname,moon_name)
    # end
    println("           Fitted posterior params from ",mcfile)
    for i=1:length(param)-1
        println(pname[i]," : ",mean(vec(par_mcmc[:,iburn:nsteps,i]))," ± ",std(vec(par_mcmc[:,iburn:nsteps,i]))) # writedlm(results,pbest_global)
    end
        println("σ_sys2 : ",mean(vec(par_mcmc[:,iburn:nsteps,end]))," ± ",std(vec(par_mcmc[:,iburn:nsteps,end])))
        println("Derived Parameters")
    for i=1:length(param)
        if i%5 == 0
            println(pname[i-4]," * M_star: ",mean(vec(par_mcmc[:,iburn:nsteps,i-4])).* CGS.MSUN/CGS.MEARTH," ± ",std(vec(par_mcmc[:,iburn:nsteps,i-4])).* CGS.MSUN/CGS.MEARTH)
            println("ecc : ",mean(sqrt.(vec(par_mcmc[:,iburn:nsteps,i]).^2 .+ vec(par_mcmc[:,iburn:nsteps,i-1]).^2))," ± ",std(sqrt.(vec(par_mcmc[:,iburn:nsteps,i]).^2 .+ vec(par_mcmc[:,iburn:nsteps,i-1]).^2)))
            println(pname[i-3]," * days: ",mean(vec(par_mcmc[:,iburn:nsteps,i-3]))," ± ",std(vec(par_mcmc[:,iburn:nsteps,i-3])))
        end     
    end
    # println(" lprob: ",mean(vec(lprob_mcmc[iburn:nsteps]))," ± ",std(vec(lprob_mcmc[iburn:nsteps])))
    sigtot=[sqrt((sqrt(mean(vec(par_mcmc[:,iburn:nsteps,end]))).*3600*24)^2 + sigma^2) ]
    println("σ_tot : ",sigtot)
end
# Print results from likelihood fit
function fit_vals(sigma::Real,nyear::Real,grid_type_nplanet::String,case_num::Int,include_moon::Bool=false)

    if case_num==1 && isfile(string("FITS/fromEMB/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2"))
        fitfile=string("FITS/fromEMB/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
        mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    elseif case_num==2 && isfile(string("FITS/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2"))
        fitfile=string("FITS/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
        # mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    else
        return println("FITS file for case ",case_num," with ",grid_type_nplanet," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
    end
    # jldmc=jldopen(mcfile,"r")
    jldfit=jldopen(String(fitfile),"r")
    pname=["mu_1","P_1","t01","ecos1","esin1",
          "mu_2","P_2","t02","ecos2","esin2",
          "mu_3","P_3","t03","ecos3","esin3", 
            "tcosϕ","tsinϕ","Δϕ","σ_sys2"]
 
    if grid_type_nplanet=="p4" || grid_type_nplanet=="p3moonp4"
        param=jldfit["best_p4"]
        lprob=jldfit["lprob_best_p4"]
        pname=[pname[1:15];"mu_4";"P_4";"t04";"ecos4";"esin4";pname[end]]
    elseif grid_type_nplanet=="p3"
        param=jldfit["best_p3"]
        lprob=jldfit["lprob_best_p3"]
		elseif grid_type_nplanet=="p2"
        param=jldfit["best_p2"]
        lprob=jldfit["lprob_best_p2"]
    end
    # if sim=="noEMB" && model=="moon"
    #         param=jldfit["best_p3"]
    #         lprob=jldfit["lprob_best_p3"]
    if include_moon && grid_type_nplanet=="p3moon"
            param=jldfit["best_dp"]
            lprob=jldfit["lprob_best_dp"]
    end
    tt0,tt,ttmodel,sigtt=jldfit["tt0"],jldfit["tt"],jldfit["ttmodel"],jldfit["sigtt"]
    nplanet,ntrans=jldfit["nplanet"],jldfit["ntrans"]
    nt1,nt2=jldfit["ntrans"][1],jldfit["ntrans"][2]
    Nobs=sum([nt1,nt2])
    jmax=5
    weight=ones(nt1+nt2)./ sigtt.^2 
    # println(param)
    chi2=0
    # Perform fit with best params, and calculate parameters covariances
    if include_moon
    fit=curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,false),tt0,tt,weight,param)
    chi2=chisquare(tt0,nplanet,ntrans,param,tt,sigtt,jmax,false)
    else
    fit=curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,true),tt0,tt,weight,param)
    chi2=chisquare(tt0,nplanet,ntrans,param,tt,sigtt,jmax,true)
    end
    cov=estimate_covar(fit)
    err=[sqrt(cov[i,j]) for i in 1:length(param), j in 1:length(param) if i==j ]
    println("           Fitted params from ",fitfile)
    #for i=1:length(param)
    #    println(pname[i]," : ",param[i]," ± ",err[i])
    #end
    println(" lprob: ",lprob)
    println(" chi: ",chi2)
    println("           Derived Parameters")
    masses = [param[i-4].* CGS.MSUN/CGS.MEARTH for i in 1:length(param) if i%5==0]
    mass_errs=[err[i-4] .* CGS.MSUN/CGS.MEARTH for i in 1:length(param) if i%5==0]
    periods = [param[i-3] for i in 1:length(param) if i%5==0]
    per_errs = [param[i-3] for i in 1:length(param) if i%5==0]
    ecc = [sqrt.(param[i].^2 .+ param[i-1].^2) for i in 1:length(param) if i%5==0]
    ecc_errs=[sqrt.(((param[i].^2 .* err[i].^2)./(param[i].^2 .+ param[i-1].^2)) .+ ((param[i].^2 * err[i-1].^2)./(param[i].^2 .+ param[i-1].^2))) for i in 1:length(param) if i%5==0 ]
    println("M_p[M⊕]=",masses," +/- ",mass_errs)
    println("Per [d]=",periods," +/- ",per_errs)
    println("eccen. =",ecc," +/- ",ecc_errs)
end
# Retrieve MCMC results after burn-in
function mc_vals(sigma::Real,nyear::Real,grid_type_nplanet::String,case_num=Int)
    # mcfile=string("MCMC/",model,"_fit",sigma,"s",nyear,"yrs.jld2")
    if case_num==1 && isfile(string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
        mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    elseif case_num==2 && isfile(string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
        mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    else
        return println("MCMC file for case ",case_num," with ",grid_type_nplanet," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
    end
    jldmc=jldopen(String(mcfile),"r")
    nwalkers,nsteps=jldmc["nwalkers"],jldmc["nsteps"]
    iburn,samples=jldmc["iburn"], jldmc["indepsamples"]
    par_mcmc=jldmc["par_mcmc"]; lprob_mcmc=jldmc["lprob_mcmc"]  ; param=jldmc["param"]
    println("           Posterior Parameters from ",mcfile)
    masses=[mean(vec(jldmc["par_mcmc"][:,iburn:end,i-4])) for i in 1:length(param) if i%5==0] .*CGS.MSUN/CGS.MEARTH
    mass_errs=[std(vec(jldmc["par_mcmc"][:,iburn:end,i-4])) for i in 1:length(param) if i%5==0] .*CGS.MSUN/CGS.MEARTH
    periods=[mean(vec(jldmc["par_mcmc"][:,iburn:end,i-3])) for i in 1:length(param) if i%5==0] 
    per_errs=[std(vec(jldmc["par_mcmc"][:,iburn:end,i-3])) for i in 1:length(param) if i%5==0] 
    ecc=[mean(sqrt.(vec(jldmc["par_mcmc"][:,iburn:nsteps,(i-1)]).^2 .+ vec(jldmc["par_mcmc"][:,iburn:nsteps,(i)]).^2)) for i in 1:length(param) if i%5==0]

		ecc_err1 = calc_quad_errs(mean(jldmc["par_mcmc"][:,iburn:end,4]),std(jldmc["par_mcmc"][:,iburn:end,4]),mean(jldmc["par_mcmc"][:,iburn:end,5]),std(jldmc["par_mcmc"][:,iburn:end,5]))
		ecc_err2 = calc_quad_errs(mean(jldmc["par_mcmc"][:,iburn:end,9]),std(jldmc["par_mcmc"][:,iburn:end,9]),mean(jldmc["par_mcmc"][:,iburn:end,10]),std(jldmc["par_mcmc"][:,iburn:end,10]))
		ecc_errs = [ecc_err1,ecc_err2]

    if grid_type_nplanet=="p3" || grid_type_nplanet=="p3moon"
			ecc_err3=calc_quad_errs(mean(jldmc["par_mcmc"][:,iburn:end,14]),std(jldmc["par_mcmc"][:,iburn:end,14]),mean(jldmc["par_mcmc"][:,iburn:end,15]),std(jldmc["par_mcmc"][:,iburn:end,15]))
			ecc_errs=[ecc_errs;ecc_err3]
    elseif grid_type_nplanet=="p4" || grid_type_nplanet=="p3moonp4"
			ecc_err4=calc_quad_errs(mean(jldmc["par_mcmc"][:,iburn:end,14]),std(jldmc["par_mcmc"][:,iburn:end,14]),mean(jldmc["par_mcmc"][:,iburn:end,15]),std(jldmc["par_mcmc"][:,iburn:end,15]))
			ecc_err3=calc_quad_errs(mean(jldmc["par_mcmc"][:,iburn:end,19]),std(jldmc["par_mcmc"][:,iburn:end,19]),mean(jldmc["par_mcmc"][:,iburn:end,20]),std(jldmc["par_mcmc"][:,iburn:end,20]))
			ecc_errs=[ecc_errs;ecc_err4;ecc_err3]
    end
		#if grid_type_nplanet=="p3moon" || grid_type_nplanet=="p3moonp4"
			#tmax_errs=calc_quad_errs(mean(jldmc["par_mcmc"][:,iburn:end,19]),std(jldmc["par_mcmc"][:,iburn:end,19]),mean(jldmc["par_mcmc"][:,iburn:end,20]),std(jldmc["par_mcmc"][:,iburn:end,20]))
		#end
    #ecc_errs= [sqrt((vec(jldmc["par_mcmc"][:,iburn:nsteps,(i-1)]).^2 .* (std(vec((jldmc["par_mcmc"][:,iburn:nsteps,(i-1)])).^2))  / (vec(jldmc["par_mcmc"][:,iburn:nsteps,(i-1)]).^2 .+ vec(jldmc["par_mcmc"][:,iburn:nsteps,(i)]).^2)) .+ (vec(jldmc["par_mcmc"][:,iburn:nsteps,(i)]).^2 .* (std(vec((jldmc["par_mcmc"][:,iburn:nsteps,(i)])).^2))  / (vec(jldmc["par_mcmc"][:,iburn:nsteps,(i-1)]).^2 .+ vec(jldmc["par_mcmc"][:,iburn:nsteps,(i)]).^2))) for i in 1:length(param) if i%5==0]

    # masses=[mean(vec(par_mcmc[:,iburn:nsteps,i-4])).*CGS.MSUN/CGS.MEARTH for i in 1:length(param) if i%5==0]
    # mass_errs=[std(vec(par_mcmc[:,iburn:nsteps,i-4])).*CGS.MSUN/CGS.MEARTH for i in 1:length(param) if i%5==0]
    # periods=[mean(vec(par_mcmc[:,iburn:nsteps,i-3])) for i in 1:length(param) if i%5==0]
    # per_errs=[std(vec(par_mcmc[:,iburn:nsteps,i-3])) for i in 1:length(param) if i%5==0]

    
    sigsys=(mean(vec(jldmc["par_mcmc"][:,iburn:end,end]))).* 3600*24
		sigsys_err=(std(vec(jldmc["par_mcmc"][:,iburn:end,end]))).* 3600*24
    sigtot=sqrt(sigsys^2 + sigma^2) 

    println("Retrieved values.")
    println("M_p[M⊕]=",masses," +/- ",mass_errs)
    println("eccen. =",ecc," +/- ",ecc_errs)
    println("Per [d]=",periods," +/- ",per_errs)
    println("σsys[s]=",sigsys," +/- ",sigsys_err)
    println("σtot[s]=",sigtot)
    # return masses, mass_errs, periods, per_errs, sigtot
end
function test_print()
    model="p3moonp4" ; case_num=2
    sigma=30 ;   nyear=30
    include_moon=false
    #fit_vals(sigma,nyear,model,case_num,include_moon)
    mc_vals(sigma,nyear,model,case_num)
    # fitfile=string("FITS/fromEMB/",grid_type_nplanet,"_fit",sigma,"s",nyear,"yrs.jld2")
    # mcfile=string("MCMC/fromEMB/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
end
