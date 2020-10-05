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
    r.MouseID != "RJ58" && # blind
    r.MouseID != "CD09" && # biting, see B1_CD09_2020-07-13 minute30
    r.MouseID != "RJ27" && # water leak
    r.MouseID != "RJ35" && # water leak
    r.MouseID != "RJ43" && # water leak
    r.MouseID != "RJ57" && # biting, see B1_RJ57_2020-09-28 minute 20:38
    !(r.MouseID in sixty_days_old) &&
    r.ProtocolSession == 1
    ,df)
end
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
cas_df = filter(r->
    r.Gen == "Rbp4-cre"
    ,Cas_s)
age_df = filter(r->
    r.Sex != "c"
    ,Age_s)
## AfterLast with Juveniles
age_afterlast = DoubleAnalysis(age_df,:Age,:AfterLast, yspan = (0,6))
age_afterlast.JarqueBera
age_afterlast.nonparametric_plot
fm1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + (1|MouseID)),age_df))
fm2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + (1|MouseID)),age_df))
fm3 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + Sex + (1|MouseID)),age_df))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
open_html_table(age_afterlast.mice_summary)
#=
Likelihood ratio test on
AfterLast ~ 1 + Age + (1|MouseID) versus AfterLast ~ 1 + (1|MouseID): p < 0.01,
effect of Juveniles in number of attempts after last reward = -0.70 ± 0.19
=#
## Correct with Juveniles
age_correct = DoubleAnalysis(age_df,:Age,:CorrectLeave, yspan = (0,1))
age_correct.JarqueBera
age_correct.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(CorrectLeave ~ 1 + (1|MouseID)),age_df))
fm2 = fit!(LinearMixedModel(@formula(CorrectLeave ~ 1 + Age + (1|MouseID)),age_df))
fm3 = fit!(LinearMixedModel(@formula(CorrectLeave ~ 1 + Age + Sex + (1|MouseID)),age_df))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
#=
Likelihood ratio test on
CorrectLeave ~ 1 + Age + (1|MouseID) versus CorrectLeave ~ 1 + (1|MouseID): p < 0.01,
effect of Juveniles on probability of correct leave = -0.07 ± 0.02
=#
## AfterLast with Caspase
cas_afterlast = DoubleAnalysis(cas_df,:Virus,:AfterLast, yspan = (0,6))
cas_afterlast.JarqueBera
cas_afterlast.nonparametric_plot
fm1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + (1|MouseID)),cas_df))
fm2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Virus + (1|MouseID)),cas_df))
Likelyhood_Ratio_test(fm1,fm2)
#=
Likelihood ratio test on
AfterLast ~ 1 + Virus + (1|MouseID) versus AfterLast ~ 1 + (1|MouseID): p < 0.01,
effect of Caspase in number of attempts after last reward = -1.26 ± 0.39
=#
## Correct with Caspase
cas_correct = DoubleAnalysis(cas_df,:Virus,:CorrectLeave, yspan = (0,1))
cas_correct.JarqueBera
cas_correct.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(CorrectLeave ~ 1 + (1|MouseID)),cas_df))
fm2 = fit!(LinearMixedModel(@formula(CorrectLeave ~ 1 + Virus + (1|MouseID)),cas_df))
Likelyhood_Ratio_test(fm1,fm2)
#=
Likelihood ratio test on
CorrectLeave ~ 1 + Age + (1|MouseID) versus CorrectLeave ~ 1 + (1|MouseID): p < 0.01,
effect of Juveniles on probability of correct leave = -0.12 ± 0.03
=#
## Interpoke with Juveniles
limit = quantile(collect(skipmissing(Age_p.PreInterpoke)),0.95)
f_age = filter(r ->
    r.PreInterpoke < limit
    ,Age_p)
age_interpoke = DoubleAnalysis(f_age,:Age,:PreInterpoke)
age_interpoke.JarqueBera
age_interpoke.nonparametric_plot
fm1 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + (1|MouseID)),f_age))
fm2 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Age + (1|MouseID)),f_age))
fm3 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Age + Sex + (1|MouseID)),f_age))
Likelyhood_Ratio_test(fm1,fm2)
#=
Likelihood ratio test on
Interpoke ~ 1 + Age + (1|MouseID) versus Interpoke ~ 1 + (1|MouseID): n.s.
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
Likelyhood_Ratio_test(fm1,fm2)
#=
Likelihood ratio test on
Interpoke ~ 1 + Virus + (1|MouseID) versus Interpoke ~ 1 + (1|MouseID): n.s.
=#
## Travel_to with Juveniles
limit = quantile(collect(skipmissing(Age_s.Travel_to)),0.95)
f_age = filter(r ->
    r.Travel_to < limit
    ,Age_s)
age_interpoke = DoubleAnalysis(f_age,:Age,:Travel_to)
age_interpoke.JarqueBera
age_interpoke.nonparametric_plot
fm1 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + (1|MouseID)),f_age))
fm2 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Age + (1|MouseID)),f_age))
fm3 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Age + Sex + (1|MouseID)),f_age))
Likelyhood_Ratio_test(fm1,fm2)
#=
Likelihood ratio test on
Travel_to ~ 1 + Age + (1|MouseID) versus Travel_to ~ 1 + (1|MouseID): n.s.
=#
## Travel_to with Caspase
limit = quantile(collect(skipmissing(Cas_s.Travel_to)),0.95)
f_cas = filter(r ->
    r.Travel_to < limit
    ,Cas_s)
cas_interpoke = DoubleAnalysis(f_cas,:Virus,:Travel_to)
cas_interpoke.JarqueBera
cas_interpoke.nonparametric_plot
fm1 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + (1|MouseID)),f_cas))
fm2 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Virus + (1|MouseID)),f_cas))
Likelyhood_Ratio_test(fm1,fm2)
#=
Likelihood ratio test on
Travel_to ~ 1 + Virus + (1|MouseID) versus Interpoke ~ 1 + (1|MouseID): n.s.
=#

[x for x in Cas_s.Travel_to]
V = Vector{Float64}(undef,length(Cas_s.Travel_to))






















## AfterLast with outlier
##
cas_afterlast = DoubleAnalysis(cas_df,:Virus,:AfterLast)
cas_afterlast.JarqueBera
cas_afterlast.parametric_plot
cas_afterlast.nonparametric_plot
# savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/AfterLastCas.pdf")
## Outlier analysis minute 30:04 starts biting and breaking the beam
o_df = filter(r -> r.Virus == "Caspase",cas_df)
 ###### After Last #######
check_cd9(o_df,:AfterLast)
 ###### InterPoke #######
check_cd9(o_df,:PreInterpoke)
 ###### Travel to #######
check_cd9(o_df,:Travel_to)
##
cd9_s = filter(r -> r.MouseID == "CD09",cas_df)

@df cd9_s scatter(:AfterLast,:PreInterpoke, markersize = 3)
@df cd9_s scatter(:AfterLast,:Travel_to, markersize = 3)

cd9_p = filter(r -> r.MouseID == "CD09",Cas_p)
@df cd9_p scatter(:PokeInStreak,:PreInterpoke, markersize = 3)
open_html_table(cd9_p)
open_html_table(cd9_s)
##
check_cd9(Cas_p,:PokeInStreak,:PreInterpoke)
##
using GLM

nullmodel = lm(@formula(AfterLast ~ 1), age_df);
model = lm(@formula(AfterLast ~ 1 + Age), age_df);
bigmodel = lm(@formula(AfterLast ~ 1 + Age + Sex), age_df);
fff = ftest(nullmodel.model, model.model,bigmodel.model)
##
age_afterlast = DoubleAnalysis(age_df,:Age,:AfterLast)
m_age = age_afterlast
##
using MixedModels
