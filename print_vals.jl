using PyPlot,JLD2,Statistics,LsqFit
import Main.TTVFaster.ttv_wrapper
# include("CGS.jl")
# Print results from MCMC run after burn-in
function print_vals(sigma::Float64,nyear::Float64,sim::String,model::String)
    mcfile=string("MCMC/",model,"_fit",sigma,"s",nyear,"yrs.jld2")
    if String(sim)=="EMB" && isfile(string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
        mcfile=string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    elseif isfile(string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
        mcfile=string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    else
        return println("MCMC file for ",sim," with ",model," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
    end
    jldmc=jldopen(String(mcfile),"r")
    par_mcmc,lprob_mcmc=jldmc["par_mcmc"],jldmc["lprob_mcmc"]
    iburn,samples=jldmc["iburn"], jldmc["indepsamples"]
    nwalkers,nsteps=jldmc["nwalkers"],jldmc["nsteps"]
    param=jldmc["param"]
    pname=["mu_1","P_1","t01","ecos1","esin1",
          "mu_2","P_2","t02","ecos2","esin2",
          "mu_3","P_3","t03","ecos3","esin3", 
            "tcosϕ","tsinϕ","Δϕ","σ_sys2"]
    if model=="p4"
        pname=[pname[1:15];"mu_4";"P_4";"t04";"ecos4";"esin4";pname[end]]
    end
    println("           Fitted posterior params from ",mcfile)
    for i=1:length(param)
        println(pname[i]," : ",mean(vec(par_mcmc[:,iburn:nsteps,i]))," ± ",std(vec(par_mcmc[:,iburn:nsteps,i]))) # writedlm(results,pbest_global)
    end
        println(pname[end]," : ",mean(vec(par_mcmc[:,iburn:nsteps,end]))," ± ",std(vec(par_mcmc[:,iburn:nsteps,end])))
        println("Derived Parameters")
    for i=1:length(param)
        if i%5 == 0
            println(pname[i-4]," * M_star: ",mean(vec(par_mcmc[:,iburn:nsteps,i-4])).* CGS.MSUN/CGS.MEARTH," ± ",std(vec(par_mcmc[:,iburn:nsteps,i-4])).* CGS.MSUN/CGS.MEARTH)
            println("ecc : ",mean(sqrt.(vec(par_mcmc[:,iburn:nsteps,i]).^2 .+ vec(par_mcmc[:,iburn:nsteps,i-1]).^2))," ± ",std(sqrt.(vec(par_mcmc[:,iburn:nsteps,i]).^2 .+ vec(par_mcmc[:,iburn:nsteps,i-1]).^2)))
        end     
    end
    println(" lprob: ",mean(vec(lprob_mcmc[iburn:nsteps]))," ± ",std(vec(lprob_mcmc[iburn:nsteps])))
end
# Print results from likelihood fit
function print_fits(sigma::Float64,nyear::Float64,sim::String,model::String,no_moon::Bool=true)
    fitfile=string("FITS/",model,"_fit",sigma,"s",nyear,"yrs.jld2")
    if String(sim)=="EMB" && isfile(string("FITS/fromEMB/",model,"_fit",sigma,"s",nyear,"yrs.jld2"))
        fitfile=string("FITS/fromEMB/",model,"_fit",sigma,"s",nyear,"yrs.jld2")
    elseif isfile(string("FITS/",model,"_fit",sigma,"s",nyear,"yrs.jld2"))
        fitfile=fitfile
    else
        return println("FITS file for ",sim," with ",model," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
    end
    jldfit=jldopen(String(fitfile),"r")
    param=jldfit["best_p3"]
    lprob=jldfit["lprob_best_p3"]
    pname=["mu_1","P_1","t01","ecos1","esin1",
          "mu_2","P_2","t02","ecos2","esin2",
          "mu_3","P_3","t03","ecos3","esin3", 
            "tcosϕ","tsinϕ","Δϕ","σ_sys2"]
    if model=="p4"
        param=jldfit["best_p4"]
        lprob=jldfit["lprob_best_p4"]
        pname=[pname[1:15];"mu_4";"P_4";"t04";"ecos4";"esin4";pname[end]]
    end
    if sim=="noEMB" && model=="moon"
            param=jldfit["best_p3"]
            lprob=jldfit["lprob_best_p3"]
    elseif model=="moon" 
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
    # Perform fit with best params, and calculate parameters covariances
    fit=curve_fit((tt0,params) -> ttv_wrapper(tt0,nplanet,ntrans,params,jmax,no_moon),tt0,tt,weight,param)
    cov=estimate_covar(fit)
    err=[sqrt(cov[i,j]) for i in 1:length(param), j in 1:length(param) if i==j ]
    println("           Fitted params from ",fitfile)
    for i=1:length(param)
        println(pname[i]," : ",param[i]," ± ",err[i])
    end
    println("Derived Parameters")
    for i=1:length(param)
        if i%5 == 0
            println(pname[i-4]," * M_star: ",param[i-4].* CGS.MSUN/CGS.MEARTH," ± ",err[i-4] .* CGS.MSUN/CGS.MEARTH)
            println("ecc : ",sqrt.(param[i].^2 .+ param[i-1].^2)," ± ",sqrt.(err[i].^2 .+ err[i-1].^2))
        end
    end
    println(" lprob: ",lprob)
end
# Retrieve MCMC results after burn-in
function get_vals(sigma::Float64,nyear::Float64,sim::String,model::String)
    mcfile=string("MCMC/",model,"_fit",sigma,"s",nyear,"yrs.jld2")
    if String(sim)=="EMB" && isfile(string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
        mcfile=string("MCMC/fromEMB/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    elseif isfile(string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2"))
        mcfile=string("MCMC/",model,"_mcmc",sigma,"s",nyear,"yrs.jld2")
    else
        return println("MCMC file for ",sim," with ",model," model at ",sigma," secs and ",nyear," yrs doesn't exist!!!!")
    end
    jldmc=jldopen(String(mcfile),"r")
    par_mcmc,lprob_mcmc=jldmc["par_mcmc"],jldmc["lprob_mcmc"]
    iburn,samples=jldmc["iburn"], jldmc["indepsamples"]
    nwalkers,nsteps=jldmc["nwalkers"],jldmc["nsteps"]
    param=jldmc["param"]
    # Build arrays of params
    masses=[mean(vec(par_mcmc[:,iburn:nsteps,i-4])).*CGS.MSUN/CGS.MEARTH for i in 1:length(param) if i%5==0]
    mass_errs=[std(vec(par_mcmc[:,iburn:nsteps,i-4])).*CGS.MSUN/CGS.MEARTH for i in 1:length(param) if i%5==0]
    periods=[mean(vec(par_mcmc[:,iburn:nsteps,i-3])) for i in 1:length(param) if i%5==0]
    per_errs=[std(vec(par_mcmc[:,iburn:nsteps,i-3])) for i in 1:length(param) if i%5==0]
    ecc=[mean(sqrt.(vec(par_mcmc[:,iburn:nsteps,i]).^2 .+ vec(par_mcmc[:,iburn:nsteps,i-1]).^2)) for i in 1:length(param) if i%5==0]
    ecc_errs=std(sqrt.(vec(par_mcmc[:,iburn:nsteps,i]).^2 .+ vec(par_mcmc[:,iburn:nsteps,i-1]).^2))) for i in 1:length(param) if i%5==0]
    sigtot=[sqrt((sqrt(mean(vec(par_mcmc[:,iburn:nsteps,end]))).*3600*24)^2 + sigma^2) ]
    return masses, mass_errs, periods, per_errs, sigtot
end