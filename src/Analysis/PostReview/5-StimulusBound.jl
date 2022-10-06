include(joinpath(pwd(),"src","Analysis","Filters.jl"))
## Test if interaction has more explanatory power
f_add_SB_Age = @formula(AfterLast ~ 1 + Num_Rewards+Age +
    (1|MouseID) + (Num_Rewards|MouseID));
f_mult_SB_Age = @formula(AfterLast ~ 1 + Num_Rewards*Age +
    (1|MouseID) + (Num_Rewards|MouseID));
add_SB_Age = fit(MixedModel,f_add_SB_Age,  Age_s; contrasts)
mult_SB_Age = fit(MixedModel,f_mult_SB_Age,  Age_s; contrasts)
MixedModels.likelihoodratiotest(add_SB_Age,mult_SB_Age)


f_add_SB_Cas = @formula(AfterLast ~ 1 + Num_Rewards+Virus +
    (1|MouseID) + (Num_Rewards|MouseID));
f_mult_SB_Cas = @formula(AfterLast ~ 1 + Num_Rewards*Virus +
    (1|MouseID) + (Num_Rewards|MouseID));
add_SB_Cas = fit(MixedModel,f_add_SB_Cas,  Cas_s; contrasts)
mult_SB_Cas = fit(MixedModel,f_mult_SB_Cas,  Cas_s; contrasts)
MixedModels.likelihoodratiotest(add_SB_Cas,mult_SB_Cas)
## plot gentle slope
plot_xy(Age_s,:Num_Rewards,:AfterLast; group = :Age, bin = false, legend = :top,
    ylims = (0,15), ylabel = "Consecutive failures", xlabel = "Rewards obtained")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","StimulusBound", "Age_StimBoung.pdf"))
plot_xy(Cas_s,:Num_Rewards,:AfterLast; group = :Virus, bin = false, legend = :top,
    ylims = (0,15), ylabel = "Consecutive failures", xlabel = "Rewards obtained")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","StimulusBound", "Cas_StimBoung.pdf"))
##
