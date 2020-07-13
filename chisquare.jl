include("ttv_wrapper.jl")
function chisquare(tt0, nplanet, ntrans, params, tt, sigtt, fixp3::Bool = false, p3_cur::Float64 = 0.0)
    chisq = 0.0
    # println(params, tt[1], sigtt[1])
    tt_model = ttv_wrapper(tt0, nplanet, ntrans, params, fixp3, p3_cur)
    for j=1:length(tt)
      chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
    end
    # println(nplanet)
    return chisq
end
