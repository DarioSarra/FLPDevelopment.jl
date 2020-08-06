#=
Objectiv: fit 2 Gammas mixture to reducing qqplot difference from data
steps:
1 get quantile from a MixtureModel
2 maximize loglikelihood Mixture data
func = TwiceDifferentiable(vars -> loglikelihood(MixtureModel(Gamma[Gamma(vars[1],vars[2]), Gamma(vars[3],vars[4])]), x),
                           ones(4); autodiff=:forward);
