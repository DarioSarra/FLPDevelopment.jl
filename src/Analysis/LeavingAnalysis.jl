using Revise, FLPDevelopment, BrowseTables, Random
plotlyjs(size=(600,600), tick_orientation = :out, grid = false,
    # linecolor = :black,
    markerstrokecolor = :black,
    thickness_scaling = 1,
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
    r.ProtocolSession == 1
    # r.Performance > 25 && no need because minimum is 31
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
for df in (Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Gen == "Rbp4-cre", df)
end
open_html_table(FLPDevelopment.summarydf(Age_s,Age_p))
open_html_table(FLPDevelopment.summarydf(Cas_s,Cas_p))
##
##
# Logistic regression of the leave choice Leave ~ Streak + PokeIn + PokeHierarchy + Group
# zscore full dataset all variable
# logscale poke In before model
# with predict or with the data
# bin In
# by mouse mean P of leave per bin
# mean group + sem per animal
# ########################### Filter the biting
############### Leaving Model
               ################ Virus

Cas_p.LogOut = log10.(Cas_p.Out)
Cas_s.LogDuration = log10.(Cas_s.Trial_duration)
transform!(Cas_p, [:Streak, :Out, :LogOut] .=> zscore)
transform!(Cas_s, :Streak .=> zscore)

#Compare AfterLast model
verbCas00 = @formula(Leave ~ 1 + Streak_zscore + Out_zscore + (1|MouseID));
LeaveCas00 = fit(MixedModel,verbCas00, Cas_p, Bernoulli())
verbCas0 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + (1|MouseID));
LeaveCas0 = fit(MixedModel,verbCas0, Cas_p, Bernoulli())
AlCas0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak_zscore + (1|MouseID)),Cas_s)
AIC_test(AlCas0, LeaveCas0)
AIC_test(LeaveCas00, LeaveCas0)

verbCas1 = @formula(Leave ~ 1 + Streak_zscore * Virus + LogOut_zscore * Virus +  (1|MouseID));
verbCas2 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore * Virus +  (1|MouseID));
verbCas3 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + Virus + (1|MouseID));

verbCas4 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +  (1+Streak_zscore+LogOut_zscore|MouseID));
LeaveCas4 = fit(MixedModel,verbCas4, Cas_p, Bernoulli())
verbCas5 = @formula(Leave ~ 1 + Streak_zscore * Virus + LogOut_zscore * Virus +  (1+Streak_zscore+LogOut_zscore|MouseID));
LeaveCas5 = fit(MixedModel,verbCas5, Cas_p, Bernoulli())
MixedModels.likelihoodratiotest(LeaveCas0,LeaveCas4)
MixedModels.likelihoodratiotest(LeaveCas4,LeaveCas5)


# Group Effect
LeaveCas1 = fit(MixedModel,verbCas1, Cas_p, Bernoulli())
LeaveCas2 = fit(MixedModel,verbCas2, Cas_p, Bernoulli())
LeaveCas3 = fit(MixedModel,verbCas3, Cas_p, Bernoulli())


MixedModels.likelihoodratiotest(LeaveCas0,LeaveCas1)
MixedModels.likelihoodratiotest(LeaveCas0,LeaveCas3)


Cas_p.M1_Leave = predict(LeaveCas1)
Cas_p.M2_Leave = predict(LeaveCas2)
Cas_p.M3_Leave = predict(LeaveCas3)
                          ################ Age
Age_p.LogOut = log10.(Age_p.Out)
Age_s.LogDuration = log10.(Age_s.Trial_duration)
transform!(Age_p, [:Streak, :Out, :LogOut] .=> zscore)
transform!(Age_s, :Streak .=> zscore)

#Compare AfterLast model
verbAge00 = @formula(Leave ~ 1 + Streak_zscore + Out_zscore + (1|MouseID));
LeaveAge00 = fit(MixedModel,verbAge00, Age_p, Bernoulli())
verbAge0 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + (1|MouseID));
LeaveAge0 = fit(MixedModel,verbAge0, Age_p, Bernoulli())
AlAge0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak_zscore + (1|MouseID)),Age_s)
AIC_test(AlAge0, LeaveAge0)
AIC_test(LeaveAge00, LeaveAge0)

verbAge1 = @formula(Leave ~ 1 + Streak_zscore * Age + LogOut_zscore * Age +  (1|MouseID));
verbAge2 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore * Age +  (1|MouseID));
verbAge3 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + Age + (1|MouseID));

verbAge4 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +  (1+Streak_zscore+LogOut_zscore|MouseID));
LeaveAge4 = fit(MixedModel,verbAge4, Age_p, Bernoulli())
verbAge5 = @formula(Leave ~ 1 + Streak_zscore * Age + LogOut_zscore * Age +  (1+Streak_zscore+LogOut_zscore|MouseID));
LeaveAge5 = fit(MixedModel,verbAge5, Age_p, Bernoulli())

MixedModels.likelihoodratiotest(LeaveAge0,LeaveAge4)
MixedModels.likelihoodratiotest(LeaveAge4,LeaveAge5)

# Group Effect
LeaveAge1 = fit(MixedModel,verbAge1, Age_p, Bernoulli())
LeaveAge2 = fit(MixedModel,verbAge2, Age_p, Bernoulli())
LeaveAge3 = fit(MixedModel,verbAge3, Age_p, Bernoulli())


MixedModels.likelihoodratiotest(LeaveAge0,LeaveAge1)
MixedModels.likelihoodratiotest(LeaveAge0,LeaveAge3)


Age_p.M1_Leave = predict(LeaveAge1)
Age_p.M2_Leave = predict(LeaveAge2)
Age_p.M3_Leave = predict(LeaveAge3)

## Figures
Age_p.LogOut = log10.(Age_p.Out)
Age_s.LogDuration = log10.(Age_s.Trial_duration)
LeaveModel = LeaveAge5
Age_bdf = FLPDevelopment.bootstrapdf(Age_p,LeaveAge5)
PLeave_Age, LeavexTrial_Age, Model_Age,
    HExp_Age, HCon_Age, HDiff_Age,
    MedianSurvival_Age, Survival_Age, Hazard_Age = Leave_plots(Age_p, Age_s, Age_bdf, filtering = true)
Age_top = plot(PLeave_Age, LeavexTrial_Age, Model_Age, layout = (1,3), size = (1850,600))
Age_mid = plot(HExp_Age, HCon_Age, HDiff_Age, layout = (1,3), size = (1850,600))
Age_bot = plot(MedianSurvival_Age, Survival_Age, Hazard_Age, layout = (1,3), size = (1850,600))
savefig(Age_top,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-1-Correlations.html"))
savefig(Age_mid,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-2-Heatmaps.html"))
savefig(Age_bot,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-3-Survival.html"))
Age_s.LeaveEvent = EventTime.(Age_s.LogDuration,trues(nrow(Age_s)))
Age_m1 = Survival.coxph(@formula(LeaveEvent ~ Streak), Age_s)
Age_m2 = Survival.coxph(@formula(LeaveEvent ~ Streak + Age), Age_s)
##
Cas_p.LogOut = log10.(Cas_p.Out)
Cas_s.LogDuration = log10.(Cas_s.Trial_duration)
Cas_bdf = FLPDevelopment.bootstrapdf(Cas_p,LeaveCas5)
PLeave_Cas, LeavexTrial_Cas, Model_Cas,
    HExp_Cas, HCon_Cas, HDiff_Cas,
    MedianSurvival_Cas, Survival_Cas, Hazard_Cas = Leave_plots(Cas_p, Cas_s, Cas_bdf; filtering = true)
Cas_top = plot(PLeave_Cas, LeavexTrial_Cas, Model_Cas, layout = (1,3), size = (1850,600))
Cas_mid = plot(HExp_Cas, HCon_Cas, HDiff_Cas, layout = (1,3), size = (1850,600))
Cas_bot = plot(MedianSurvival_Cas, Survival_Cas, Hazard_Cas, layout = (1,3), size = (1850,600))

savefig(Cas_top,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-1-Correlations.html"))
savefig(Cas_mid,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-2-Heatmaps.html"))
savefig(Cas_bot,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-3-Survival.html"))
Cas_s.LeaveEvent = EventTime.(Cas_s.LogDuration,trues(nrow(Cas_s)))
Cas_m1 = Survival.coxph(@formula(LeaveEvent ~ Streak), Cas_s)
Cas_m2 = Survival.coxph(@formula(LeaveEvent ~ Streak + Virus), Cas_s)
##
FAge_p, FAge_s = filter_pokestreak(Age_p)
FAge_p.LogOut = log10.(FAge_p.Out)
FAge_s.LogDuration = log10.(FAge_s.Trial_duration)
PLeave_FAge, NPokes_FAge, Model_FAge, HExp_FAge, HCon_FAge, HDiff_FAge = Leave_plots(FAge_p, FAge_s; model_plt = Model_FAge, filtering = true)
FAge_top = plot(PLeave_FAge, NPokes_FAge, Model_FAge, layout = (1,3),size = (1850,600))
FAge_bot = plot(HExp_FAge, HCon_FAge, HDiff_FAge, layout = (1,3), size = (1850,600))
savefig(FAge_top,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-4-Filt-DurCorrelations.html"))
savefig(FAge_bot,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-5-Filt-FullHeatmaps.html"))
##
##
FCas_p, FCas_s = filter_pokestreak(Cas_p)
FCas_p.LogOut = log10.(FCas_p.Out)
FCas_s.LogDuration = log10.(FCas_s.Trial_duration)
PLeave_FCas, NPokes_FCas, Model_FCas, HExp_FCas, HCon_FCas, HDiff_FCas = Leave_plots(FCas_p, FCas_s)#; model_plt = Model_FCas, filtering = true)
FCas_top = plot(PLeave_FCas, NPokes_FCas, Model_FCas, layout = (1,3), size = (1850,600))
FCas_bot = plot(HExp_FCas, HCon_FCas, HDiff_FCas, layout = (1,3), size = (1850,600))
savefig(FCas_top,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-4-Filt-DurCorrelations.html"))
savefig(FCas_bot,joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-5-Filt-FullHeatmaps.html"))
######################################## Survival Analysis
##
Age_s.LogDuration = log10.(Age_s.Trial_duration)
Age_MedianSurvival = mediansurvival_analysis(Age_s,:LogDuration, :Age)
Age_Survival = function_analysis(Age_s,:LogDuration, survivalrate_algorythm; grouping = :Age)
plot!(Age_Survival, xlabel = "Time (log10 s)", ylabel = "Survival rate", label = "")
Age_Hazard = function_analysis(Age_s,:LogDuration, hazardrate_algorythm, grouping = :Age)
plot!(Age_Hazard, xlabel = "Time (log10 s)", ylabel = "Hazard rate")
plot(Age_MedianSurvival, Age_Survival,Age_Hazard)
Age_s.LeaveEvent = EventTime.(Age_s.LogDuration,trues(nrow(Age_s)))
Age_m1 = Survival.coxph(@formula(LeaveEvent ~ Streak), Age_s)
Age_m2 = Survival.coxph(@formula(LeaveEvent ~ Streak + Age), Age_s)
##
Cas_s.LogDuration = log10.(Cas_s.Trial_duration)
Cas_MedianSurvival = mediansurvival_analysis(Cas_s,:LogDuration, :Virus)
Cas_Survival = function_analysis(Cas_s,:LogDuration, survivalrate_algorythm; grouping = :Virus)
Cas_Hazard = function_analysis(Cas_s,:LogDuration, hazardrate_algorythm, grouping = :Virus)
plot(Cas_MedianSurvival, Cas_Survival, Cas_Hazard)
Cas_s.LeaveEvent = EventTime.(Cas_s.LogDuration,trues(nrow(Cas_s)))
Cas_m1 = Survival.coxph(@formula(LeaveEvent ~ Streak), Cas_s)
Cas_m2 = Survival.coxph(@formula(LeaveEvent ~ Streak + Virus), Cas_s)
##
fit(KaplanMeier, Cas_s.LogDuration, trues(nrow(Cas_s)))
EventTime.(Cas_s.Streak)
fit(CoxModel,EventTime.(Cas_s.Streak), Cas_s.LogDuration)
)
