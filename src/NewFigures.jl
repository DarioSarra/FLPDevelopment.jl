# Afterlast over trials
##
using Revise, FLPDevelopment, BrowseTables
gr(size=(600,600), tick_orientation = :out, grid = false,
    linecolor = :black,
    markerstrokecolor = :black,
    thickness_scaling = 2,
    markersize = 8)
##
include("Young_to_run2.jl")
include("Caspase_to_run.jl")
for df in (Age_p, Age_b, Age_s, Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Protocol == "90/90" &&
    r.MouseID != "CD09" && # biting, see B1_CD09_2020-07-13 minute30
    r.MouseID != "RJ58" && # blind
    r.MouseID != "RJ67" && # biting, see B3_RJ67_2020-09-28 minute 7:33
    !(r.MouseID in first_females_group) &&
    # r.Performance > 25 &&
    r.ProtocolSession == 1
    # r.Streak < 75 && #checking
    # previously tried filters
    # r.MouseID != "RJ27" && # water leak
    # r.MouseID != "RJ35" && # water leak
    # r.MouseID != "RJ43" && # water leak
    # r.MouseID != "RJ57" && # biting, see B1_RJ57_2020-09-28 minute 20:38
    # r.MouseID != "RJ70" && # biting, see B1_RJ70_2020-09-28 minute 24:23
    # !(r.MouseID in second_females_juveniles) &&
    # !(r.MouseID in sixty_days_old) &&
    ,df)
end
##### Selection criteria ##
# AL over trials

Cas_s[!,:BinnedStreak] = bin_axis(Cas_s.Streak; unit_step = 4)
res3 = summary_xy(Cas_s,:BinnedStreak,:AfterLast; group = :Virus)
@df res3 plot(string.(:BinnedStreak),:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xrotation = 50)
# trials over time
Cas_s.Start
Cas_s[!,:BinnedStart] = bin_axis(Cas_s.Start./60; unit_step = 2)
res3 = summary_xy(Cas_s,:BinnedStart,:Streak; group = :Virus)
@df res3 plot(string.(:BinnedStart),:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, legend = :topleft)
# AL normalisation
df1 = group_kde(Cas_s,:AfterLast; group = [:Virus], points = 50)
@df df1 plot(:Xaxis,:Mean, ribbon = :Sem, xlims = (0,25), linewidth = 1, group = :Virus)
#####
##
#=
to do list
- group median and ci plus mean animal afterlast plot
- single animal distributions
- CD9 interpoke, trial duaration and travel time against Rbp4
- CD9 long trials interpoke, trial duaration and travel time against short trials
- CD9 video
=#
## Afterlast df selection
limit_cas = quantile(collect(skipmissing(Cas_s.AfterLast)),0.95)
limit_age = quantile(collect(skipmissing(Cas_s.AfterLast)),0.95)
cas_df = filter(r->
    r.Trial_duration < 30 &&
    r.AfterLast < 8 &&
    r.Gen == "Rbp4-cre"
    ,Cas_s)
age_df = filter(r->
    r.Trial_duration < 30 &&
    r.AfterLast < 8
    ,Age_s)
## AfterLast with Juveniles
age_afterlast = DoubleAnalysis(age_df,:Age,:AfterLast)
age_afterlast.JarqueBera
age_afterlast.nonparametric_plot
confint(age_afterlast.UnequalVarianceT)
age_afterlast.parametric_plot
age_afterlast.parametric_summary
res = age_afterlast.mice_summary
@df res boxplot(:Age,:AfterLast)

fm1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + (1|MouseID)),age_df))
fm2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + (1|MouseID)),age_df))
fm3 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + Sex + (1|MouseID)),age_df))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)

fm4 = fit(MixedModel,@formula(AfterLast ~ 1 + (1|MouseID)),age_df,Poisson())
fm5 = fit(MixedModel,@formula(AfterLast ~ 1 + Age + (1|MouseID)),age_df,Poisson())
fm6 = fit(MixedModel,@formula(AfterLast ~ 1 + Age + Trial_duration + (1|MouseID)),age_df,Poisson())
age_df
Likelyhood_Ratio_test(fm5,fm6)
scatter(age_df.AfterLast,predict(fm5),markersize=2, legend  = false)
p = predict(fm6)
y = age_df.AfterLast
r = @. sign(y - p) * sqrt(devresid(Poisson(), y, p))
scatter(r,age_df.AfterLast;markersize=1, stroke = 1, legend = false)
qqnorm(r; qqline = :R, markersize = 2)
fitted(fm5)
# open_html_table(age_afterlast.mice_summary)
# open_html_table(filter(r-> r.MouseID == "RJ67", age_df))
#=
Likelihood ratio test on
AfterLast ~ 1 + Age + (1|MouseID) versus AfterLast ~ 1 + (1|MouseID): p < 0.01,
effect of Juveniles in number of attempts after last reward = -0.70 ± 0.19
=#
## Incorrect with Juveniles
age_correct = DoubleAnalysis(age_df,:Age,:IncorrectLeave, yspan = (0,1))
age_correct.JarqueBera
age_correct.nonparametric_plot
age_correct.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + (1|MouseID)),age_df))
fm2 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + Age + (1|MouseID)),age_df))
fm3 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + Age + Sex + (1|MouseID)),age_df))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
age_correct.parametric_summary
res = age_correct.mice_summary
res[!,:Sex] = [x in females ? "F" : "M" for x in res.MouseID]
res[!,:Combo] = res.Sex .* res.Age
group_summary(res,:Combo, :IncorrectLeave)
JarqueBeraTest(residuals(fm2))
# open_html_table(age_correct.mice_summary)
#=
Likelihood ratio test on
IncorrectLeave ~ 1 + Age + (1|MouseID) versus IncorrectLeave ~ 1 + (1|MouseID): p < 0.01,
effect of Juveniles on probability of correct leave = -0.07 ± 0.02
=#
## Interpoke with Juveniles
limit = quantile(collect(skipmissing(Age_p.PreInterpoke)),0.95)
f_age = filter(r ->
    r.PreInterpoke < limit
    ,Age_p)
age_interpoke = DoubleAnalysis(f_age,:Age,:PreInterpoke)
age_interpoke.JarqueBera
age_interpoke.nonparametric_plot
age_interpoke.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + (1|MouseID)),f_age))
fm2 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Age + (1|MouseID)),f_age))
fm3 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Age + Sex + (1|MouseID)),f_age))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
age_interpoke.parametric_summary
#=
Likelihood ratio test on
Interpoke ~ 1 + Age + (1|MouseID) versus Interpoke ~ 1 + (1|MouseID): n.s.
=#
## Travel_to with Juveniles
limit = quantile(collect(skipmissing(Age_s.Travel_to)),0.95)
limit = median(Age_s.Travel_to)
f_age = filter(r ->
    r.Travel_to < limit
    ,Age_s)
age_travel = DoubleAnalysis(f_age,:Age,:Travel_to)
age_travel.JarqueBera
age_travel.nonparametric_plot
age_travel.parametric_plot
f_age[!,:Travel_to] = map(Float64,f_age.Travel_to)
fm1 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + (1|MouseID)),f_age))
fm2 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Age + (1|MouseID)),f_age))
fm3 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Age + Sex + (1|MouseID)),f_age))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
age_travel.parametric_summary
#=
Likelihood ratio test on
Travel_to ~ 1 + Age + (1|MouseID) versus Travel_to ~ 1 + (1|MouseID): n.s.
=#
## AfterLast with Caspase
cas_afterlast = DoubleAnalysis(cas_df,:Virus,:AfterLast, yspan = (0,6))
cas_afterlast.JarqueBera
cas_afterlast.nonparametric_plot
cas_afterlast.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + (1|MouseID)),cas_df))
fm2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Virus + (1|MouseID)),cas_df))
Likelyhood_Ratio_test(fm1,fm2)
cas_afterlast.mice_summary
cas_afterlast.parametric_summary
#=
Likelihood ratio test on
AfterLast ~ 1 + Virus + (1|MouseID) versus AfterLast ~ 1 + (1|MouseID): p < 0.01,
effect of Caspase in number of attempts after last reward = -1.26 ± 0.39
=#
## Incorrect with Caspase
cas_correct = DoubleAnalysis(cas_df,:Virus,:IncorrectLeave, yspan = (0,1))
cas_correct.JarqueBera
cas_correct.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + (1|MouseID)),cas_df))
fm2 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + Virus + (1|MouseID)),cas_df))
Likelyhood_Ratio_test(fm1,fm2)
cas_correct.parametric_summary
#=
Likelihood ratio test on
IncorrectLeave ~ 1 + Age + (1|MouseID) versus IncorrectLeave ~ 1 + (1|MouseID): p < 0.01,
effect of Juveniles on probability of correct leave = -0.12 ± 0.03
=#

## Interpoke with Caspase
limit = quantile(collect(skipmissing(Cas_p.PreInterpoke)),0.95)
f_cas = filter(r ->
    r.PreInterpoke < limit
    ,Cas_p)
cas_interpoke = DoubleAnalysis(f_cas,:Virus,:PreInterpoke)
cas_interpoke.JarqueBera
cas_interpoke.nonparametric_plot
fm1 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + (1|MouseID)),f_cas))
fm2 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Virus + (1|MouseID)),f_cas))
fm3 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Virus + Sex + (1|MouseID)),f_cas))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
cas_interpoke.parametric_summary
#=
Likelihood ratio test on
Interpoke ~ 1 + Virus + (1|MouseID) versus Interpoke ~ 1 + (1|MouseID): n.s.
=#

## Travel_to with Caspase
@df f_cas density(:Travel_to)
limit = quantile(collect(skipmissing(Cas_s.Travel_to)),0.95)
limit = median(Cas_s.Travel_to)
f_cas = filter(r ->
    r.Travel_to < limit
    ,Cas_s)
cas_travel = DoubleAnalysis(f_cas,:Virus,:Travel_to)
cas_travel.JarqueBera
cas_travel.nonparametric_plot
cas_travel.parametric_plot
f_cas[!,:Travel_to] = map(Float64,f_cas.Travel_to)
fm1 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + (1|MouseID)),f_cas))
fm2 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Virus + (1|MouseID)),f_cas))
fm3 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Virus + Sex + (1|MouseID)),f_cas))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
cas_travel.parametric_summary
#=
Likelihood ratio test on
Travel_to ~ 1 + Virus + (1|MouseID) versus Interpoke ~ 1 + (1|MouseID): n.s.
=#
## Number of Reward
age_rew = DoubleAnalysis(age_df,:Age,:Num_Rewards, summary_opt = :SUM)
age_rew.JarqueBera
age_rew.nonparametric_plot
age_rew.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + (1|MouseID)),age_df))
fm2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + (1|MouseID)),age_df))
fm3 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + Sex + (1|MouseID)),age_df))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
age_rew.parametric_summary
res = age_rew.mice_summary
res[!,:Sex] = [x in females ? "F" : "M" for x in res.MouseID]
res[!,:Combo] = res.Sex .* res.Age
group_summary(res,:Combo, :AfterLast)
open_html_table(age_df)
