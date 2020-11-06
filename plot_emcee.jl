using PyPlot

function plot_emcee(include_moon::Bool=false)
    figsize=(9,5)
    for i=1:5
    subplot(5,1,i)
    for j=1:nwalkers 
    plot(par_mcmc[j,1:nsteps,i])
    ylabel(pname[i])
    end
    tight_layout()
    end
    name = string("IMAGES/MCMCstepsp1.png")
    # savefig(name)
    clf()
    figsize=(9,5)
    for i=1:5
    subplot(5,1,i)
    for j=1:nwalkers 
    plot(par_mcmc[j,1:nsteps,i+5])
    ylabel(pname[i+5])
    end
    tight_layout()
    end
    name = string("IMAGES/MCMCstepsp2.png")
    # savefig(name)
    clf()
    figsize=(9,5)
    for i=1:5
    subplot(5,1,i)
    for j=1:nwalkers 
    plot(par_mcmc[j,1:nsteps,i+10])
    ylabel(pname[i+10])
    end
    tight_layout()
    end
    name = string("IMAGES/MCMCstepsp3.png")
    # savefig(name)
    clf()
    if include_moon
        figsize=(5,3)
        for i=1:3
        subplot(3,1,i)
        for j=1:nwalkers 
        plot(par_mcmc[j,1:nsteps,i+15])
        ylabel(pname[i+15])
        end
        # tight_layout()
        end
        name = string("IMAGES/MCMCstepsmoon.png")
        # savefig(name)
    end
end
# Make plot of MCMC parameters after burn-in
function plot_parameters(xvalue,yvalue)

# figsize=(8,6)
# for i=2:nparam
#   for j=1:i-1
#     scatter(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]),vec(par_mcmc[1:nwalkers,iburn:nsteps,j]))
#     xlabel(pname[i])
#     ylabel(pname[j])
#   end
# end
# name = string("IMAGES/MCMCparams",label,".png")
end