include(joinpath(pwd(),"src","Analysis","Filters.jl"))
##
NumPokes_f = @formula(Leave ~ 1 + PokeInStreak +  (1|MouseID) + (PokeInStreak|MouseID));
NumPokes_m = fit(MixedModel,NumPokes_f, Age_p, Bernoulli();contrasts)
PokeTime_f = @formula(Leave ~ 1 + LogOut_zscore +  (1|MouseID) + (LogOut_zscore|MouseID));
PokeTime_m = fit(MixedModel,PokeTime_f, Age_p, Bernoulli();contrasts)
AIC_test(PokeTime_m, NumPokes_m)
aic(PokeTime_m)
aic(NumPokes_m)
