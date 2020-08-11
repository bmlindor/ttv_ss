# plot all relevant data 
using PyPlot, JLD2
include("decompose_ttvs.jl")
label = "test"
@load "OUTPUTS/p3_fittestparams.jld2" #param_p3 lprob_p3 lprob_best pbest ntrans nplanet tt0 tt ttmodel sigtt p3in p3out np3 nphase
#@load "OUTPUTS" #mcmc 

pair_ttvs = decompose_ttvs(nplanet, ntrans, pbest)

time1 = collect(pbest[3] .+ range(0,stop=ntrans[1]-1,length=ntrans[1]) .* pbest[2])
time2 = collect(pbest[8] .+ range(0,stop=ntrans[2]-1,length=ntrans[2]) .* pbest[7])

# fig, axes =  subplots(1,2)
# ax1 = axes[1]
subplot(211)
plot(ttmodel[1:ntrans[1]],pair_ttvs[1,2,1:ntrans[1]],label=)
plot(ttmodel[1:ntrans[1]],pair_ttvs[1,3,1:ntrans[1]],label=)
errorbar(ttmodel[1:ntrans[1]],tt[1:ntrans[1]].-time1, sigtt[1:ntrans[1]],fmt='.')
# ax1.set_xlabel("Venus")
# ax2 = axes[2]
subplot(212)
plot(ttmodel[ntrans[1]+1:ntrans[1]+ntrans[2]],pair_ttvs[2,1,1:ntrans[2]],label="")
plot(ttmodel[ntrans[1]+1:ntrans[1]+ntrans[2]],pair_ttvs[2,3,1:ntrans[2]],label="")
errorbar(ttmodel[ntrans[1]+1:ntrans[1]+ntrans[2]],tt[ntrans[1]+1:ntrans[1]+ntrans[2]].-time2, sigtt[ntrans[1]+1:ntrans[1]+ntrans[2]],fmt='.')
xlabel("Earth")
legend()


