##### EV runs
#julia full_run.jl Hpp mcmc fromEV 20000 > OUTPUTS/mcmc_p2_EV.out &
#julia full_run.jl Hppp mcmc fromEV 30000 > OUTPUTS/mcmc_p3_EV.out &
julia full_run.jl Hppmp mcmc fromEV 40000 > OUTPUTS/mcmc_p3m_EV.out &
julia full_run.jl Hppmpp mcmc fromEV 50000 > OUTPUTS/mcmc_p4m_EV.out &
##### EMB runs
#julia full_run.jl Hpp mcmc fromEMB 10000 > OUTPUTS/mcmc_p2_EMB.out &
#julia full_run.jl Hppp mcmc fromEMB 20000 > OUTPUTS/mcmc_p3_EMB.out &
julia full_run.jl Hpppp mcmc fromEMB 50000 > OUTPUTS/mcmc_p4_EMB.out &
