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
############### Leaving Model
               ################ Virus

Cas_p.LogOut = log10.(Cas_p.Out)
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
                          ################ Virus
#Row data P leaving
transform!(Cas_p, :LogOut => (x -> bin_axis(x; length = 20)) => :Bin_LogOut)
df0 = filter(r -> r.Streak <= 70, Cas_p)
PLeave_Cas = FLPDevelopment.P_Leave(Cas_p,:Bin_LogOut,:Leave)
# savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-1-LogPleave.pdf"))
## PokeN over trials
Cas_s[!,:BinnedStreak] = bin_axis(Cas_s.Streak; unit_step = 4)
res1 = summary_xy(Cas_s,:BinnedStreak,:Num_pokes; group = :Virus)
NPokes_Cas = @df res1 plot(string.(:BinnedStreak),:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xlabel = "Trial", ylabel = "Number of pokes",size=(600,600))
# savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-2-PokesTrial.pdf"))
## Model Plot
rng = MersenneTwister(1234321)
Cas_samp1 = parametricbootstrap(rng,100,LeaveCas1)
Cas_sampdf = DataFrame(Cas_samp1.allpars)
Cas_bootdf = combine(groupby(Cas_sampdf,[:type, :group, :names]), :value => shortestcovint => :interval)
Cas_bootdf.coef = push!(coef(LeaveCas1), mean(ranef(LeaveCas1)[1]))
Cas_bootdf.variable = ["Intercept", "Trial", "Virus:Caspase", "Poke-time",
    "Trial & Virus: Caspase",
    "Poke-time & Virus:Caspase", "MouseID"]
transform!(Cas_bootdf, [:coef, :interval] => ByRow((c,e) -> (c -e[1], e[2]-c)) => :err)
Model_Cas = @df Cas_bootdf[1:end-1,:] scatter(:coef ,1:nrow(Cas_bootdf)-1,
    xerror = :err, xlabel = "Coefficient estimate",
    yticks = (1:nrow(Cas_bootdf), :variable))
vline!([0], linecolor = :red, legend = false)
# savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-3-LogModel.pdf"))

## P leave Heatmaps
transform!(Cas_p, :LogOut => (x -> bin_axis(x; length = 30)) => :Bin_LogOut)
transform!(Cas_p, :Streak => (x -> bin_axis(x; unit_step = 10)) => :Bin_Streak)
heat = Heatmap_group(Cas_p,:Bin_LogOut,:Bin_Streak,:Leave)
# filter!(r -> 11<= r.Bin_Streak <= 61 && 0.1 <= r.Bin_LogOut <= 1.9, heat)
heatcas = filter(r -> r.Virus == "Caspase", heat)
HExp_Cas = Heatmap_plot(heatcas,:Bin_LogOut,:Bin_Streak,:Leave)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-4-Full-HeatmapCaspaseRaw.pdf"))
##
heattom = filter(r -> r.Virus == "tdTomato", heat)
HCon_Cas = FLPDevelopment.Heatmap_plot(heattom,:Bin_LogOut,:Bin_Streak,:Leave)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-5-Full-HeatmapTomatoRaw.pdf"))
##
diff = Heatmap_difference(heat,:Bin_LogOut,:Bin_Streak,:Leave; grouping = :Virus, adjust = :trim)
HDiff_Cas = Heatmap_plot(diff,:Bin_LogOut,:Bin_Streak,:Leave; colorlims = (-1, 1))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Cas-6-Full-HeatmapDiffRaw.pdf"))
##
                          ################ Age
#Row data P leaving
transform!(Age_p, :LogOut => (x -> bin_axis(x; length = 20)) => :Bin_LogOut)
# df0 = filter(r -> r.Streak <= 70, Age_p)
PLeave_Age = FLPDevelopment.P_Leave(Age_p,:Bin_LogOut,:Leave)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-1-LogPleave.pdf"))
## PokeN over trials
Age_s[!,:BinnedStreak] = bin_axis(Age_s.Streak; unit_step = 4)
res1 = summary_xy(Age_s,:BinnedStreak,:Num_pokes; group = :Age)
NPokes_Age = @df res1 plot(string.(:BinnedStreak),:Mean, group = :Age, linecolor = :auto,
    ribbon = :Sem, xlabel = "Trial", ylabel = "Number of pokes",size=(600,600))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-2-PokesTrial.pdf"))
## Model Plot
rng = MersenneTwister(1234321)
Age_samp1 = parametricbootstrap(rng,100,LeaveAge1)
Age_sampdf = DataFrame(Age_samp1.allpars)
Age_bootdf = combine(groupby(Age_sampdf,[:type, :group, :names]), :value => shortestcovint => :interval)
Age_bootdf.coef = push!(coef(LeaveAge1), mean(ranef(LeaveAge1)[1]))
Age_bootdf.variable = ["Intercept", "Trial", "Age:Juveniles", "Poke-time",
    "Trial & Age: Juveniles",
    "Poke-time & Age:Juveniles", "MouseID"]
transform!(Age_bootdf, [:coef, :interval] => ByRow((c,e) -> (c -e[1], e[2]-c)) => :err)
Model_Age = @df Age_bootdf[1:end-1,:] scatter(:coef ,1:nrow(Age_bootdf)-1,
    xerror = :err, xlabel = "Coefficient estimate",
    yticks = (1:nrow(Age_bootdf), :variable))
vline!([0], linecolor = :red, legend = false)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-3-LogModel.pdf"))
## P leave Heatmaps
transform!(Age_p, :LogOut => (x -> bin_axis(x; length = 30)) => :Bin_LogOut)
transform!(Age_p, :Streak => (x -> bin_axis(x; unit_step = 10)) => :Bin_Streak)
Age_heat = Heatmap_group(Age_p,:Bin_LogOut,:Bin_Streak,:Leave)
filter!(r -> 11<= r.Bin_Streak <= 61 && 0.1 <= r.Bin_LogOut <= 2.1, Age_heat)
Age_heatcas = filter(r -> r.Age == "Juveniles", Age_heat)
HExp_Age = Heatmap_plot(Age_heatcas,:Bin_LogOut,:Bin_Streak,:Leave)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-4-HeatmapJuvenilesRaw.pdf"))
##
Age_heattom = filter(r -> r.Age == "Adults", Age_heat)
HCon_Age = FLPDevelopment.Heatmap_plot(Age_heattom,:Bin_LogOut,:Bin_Streak,:Leave)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-5-HeatmapAdultsRaw.pdf"))
##
Age_diff = Heatmap_difference(Age_heat,:Bin_LogOut,:Bin_Streak,:Leave; grouping = :Age, adjust = :trim)
HDiff_Age = Heatmap_plot(Age_diff,:Bin_LogOut,:Bin_Streak,:Leave; colorlims = (-1, 1))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","Age-6-HeatmapDiffRaw.pdf"))
