#=
Objectiv: fit 2 Gammas mixture to reducing qqplot difference from data
steps:
1 loglikelihood Mixture data
2 qqplot mixture distribution
=#
##
function custom_func(x,params)
    if any(params .< 0)
        return 1000000
    else
        return - loglikelihood(MixtureModel(Gamma[Gamma(params[1],params[2]), Gamma(params[3],params[4])]), x)
    end
end
# opt = optimize(vars -> custom_func(x,vars), ones(4))
# res = opt.minimizer

function mixture_gamma(x)
    opt = optimize(vars -> custom_func(x,vars), ones(4))
    res = opt.minimizer
    MixtureModel(Gamma[Gamma(res[1],res[2]), Gamma(res[3],res[4])])
end
##
function custom_func_weighted(x,params)
    if any(params .< 0)
        return 100000000
    elseif !(0 < params[5] <1)
        return 100000000
    else
        return - loglikelihood(MixtureModel(Gamma[Gamma(params[1],params[2]), Gamma(params[3],params[4])],[params[5],1 - params[5]]), x)
    end
end
# opt = optimize(vars -> custom_func(x,vars), ones(4))
# res = opt.minimizer

function mixture_gamma_weighted(x)
    opt = optimize(vars -> custom_func(x,vars), ones(5))
    res = opt.minimizer
    MixtureModel(Gamma[Gamma(res[1],res[2]), Gamma(res[3],res[4])],[res[5],1 - res[5]])
end
##
function cust_mix_expo(x,params)
    if any(params .< 0)
        return 100000000
    elseif !(0 < params[3] <1)
        return 100000000
    else
        return - loglikelihood(MixtureModel(Exponential[Exponential(params[1]), Exponential(params[2]),[params[3],1-params[3]]]), x)
    end
end

function mixture_exp_weighted(x)
    opt = optimize(vars -> cust_mix_expo(x,vars), ones(3))
    res = opt.minimizer
    MixtureModel(Exponential[Exponential(res[1]), Exponential(res[2])],[res[3],1-res[3]])
end
##
function cust_expo(x,params)
    if any(params .< 0)
        return 100000000
    else
        return - loglikelihood(MixtureModel(Gamma[Exponential(params[1]), Gamma(params[2],params[3])]), x)
        # return - loglikelihood(MixtureModel(Exponential[Exponential(params[1]), Exponential(params[2])]), x)
    end
end

function mixture_exp(x)
    opt = optimize(vars -> cust_expo(x,vars), ones(3))
    res = opt.minimizer
    MixtureModel(Gamma[Exponential(res[1]), Gamma(res[2],res[3])])
    # MixtureModel(Exponential[Exponential(res[1]), Exponential(res[2])])

end
