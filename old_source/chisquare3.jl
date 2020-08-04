include("ttv_wrapper3.jl")
function chisquare3(param)
    chisq = 0.0
    tt_model = ttv_wrapper3(tt0,param)
    for j=1:length(tt0)
      chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
    end
    return chisq
end
