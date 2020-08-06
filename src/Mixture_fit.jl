#=
Objectiv: fit 2 Gammas mixture to reducing qqplot difference from data
steps:
1 get quantile from a MixtureModel
2 maximize loglikelihood Mixture data
=#
##
x = pre_cas.Travel_to
function custom_func(params)
    if any(params .< 0)
        return 1000000
    else
        return - loglikelihood(MixtureModel(Gamma[Gamma(params[1],params[2]), Gamma(params[3],params[4])]), x)
    end
end
opt = optimize(custom_func, ones(4))
res = opt.minimizer
##
check = MixtureModel(Gamma[Gamma(res[1],res[2]), Gamma(res[3],res[4])])
check1 = check.components[1]
scaling1 = m/maximum(pdf(check1,0:60))
check2 = check.components[2]
check3 = fit(Gamma,x)
scaling3 = m/maximum(pdf(check3,0:60))
##
p = histogram(x, nbins = 100)
hist_y = p.series_list[1].plotattributes[:y]
m = maximum(filter(!isnan,hist_y))

scaling = m/maximum(pdf(check,0:60))
plot!(0:60,pdf(check,0:60)*scaling)

scaling1 = 50/maximum(pdf(check1,0:60))
plot!(0:60,pdf(check1,0:60)*scaling1, linecolor = :cyan)

scaling2 = m/maximum(pdf(check2,0:60))
plot!(0:60,pdf(check2,0:60)*scaling2, linecolor = :red)
##
p = histogram(x, nbins = 100, color = :grey, fillalpha = 0.3, linewidth = 0)
histogram!(rand(check,length(x)), nbins = 100, fillalpha = 0.3, linewidth = 0)
histogram!(rand(check1,Int64(round(length(x)/2))), nbins = 100, fillalpha = 0.3, linewidth = 0)
histogram!(rand(check2,Int64(round(length(x)/2))), nbins = 100, fillalpha = 0.3, linewidth = 0)

##
lims = (0,0.135)
plot(0:60,pdf(check,0:60),ylims = lims, linecolor = :magenta, label = "Convex combination")
plot!(0:60,pdf(check2,0:60),ylims = lims, linecolor = :cyan, xticks = 0:5:60, label = "First component")
plot!(0:60,pdf(check1,0:60),ylims = lims, linecolor = :red, label = "Second component")
q95 = quantile(check2,0.95)
vline!([q95], label = "95th percentile\nfirst component")
