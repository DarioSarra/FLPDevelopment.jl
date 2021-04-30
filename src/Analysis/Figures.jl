using Revise, FLPDevelopment, BrowseTables
##
# Font size and type
xyfont = font(18, "Times")
legfont = font(32, "Times")

# Figure/Plot image size
pixel_factor = 600
ps_x = pixel_factor * 1
ps_y = pixel_factor * 1
fig_size = (ps_x, ps_y)

theme(:default)
pgfplotsx(size = fig_size,
    tick_orientation = :out,
    grid = false,
    markerstrokecolor = :black,
    markersize = 8,
    thickness_scaling = 1.5,
    titlefont = xyfont,
    legendfont = legfont)
plot()
PGFPlotsX.latexengine()
##
include("Young_to_run2.jl")
include("Caspase_to_run.jl")
for df in (Age_p, Age_b, Age_s, Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Protocol == "90/90" &&
    r.MouseID != "CD09" && # biting, see B1_CD09_2020-07-13 minute30
    r.MouseID != "RJ58" && # blind
    r.MouseID != "RJ67" && # biting, see B3_RJ67_2020-09-28 minute 7:33
    !(r.MouseID in first_females_group) &&
    # r.Streak >= 10 && #checking
    r.ProtocolSession == 1
    # r.Performance > 25 && no need because minimum is 31
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
# fAge_p = filter(r->r.PokeDur > 0.3 &&
#     (r.Reward || ismissing(r.PostInterpoke) || (r.PostInterpoke > 0.1)) &&
#     (r.Reward || r.PreInterpoke == 0 || ismissing(r.PreInterpoke) || (r.PreInterpoke > 0.1)),
#     Age_p)
#     gd = groupby(fAge_p,[:MouseID,:Session,:Age,:Sex])
#     agedf = combine(gd) do dd
#         process_streaks(dd)
#     end
# fCas_p = filter(r->r.PokeDur > 0.3 &&
#     (r.Reward || ismissing(r.PostInterpoke) || (r.PostInterpoke > 0.1)) &&
#     (r.Reward || r.PreInterpoke == 0 || ismissing(r.PreInterpoke) || (r.PreInterpoke > 0.1)),
#     Cas_p)
#     gd = groupby(fCas_p,[:MouseID,:Session,:Virus])
#     casdf = combine(gd) do dd
#         process_streaks(dd)
#     end
# agedf[!,:BinnedStreak] = bin_axis(agedf.Streak; unit_step = 5)
# casdf[!,:BinnedStreak] = bin_axis(casdf.Streak; unit_step = 5)
#=
    Figure 1:
        - Example Session
        - Medain survival
        - Survival rate
        - Leaving time over trial
        - Model bootstrap
        - AfterLast scatter
        - Incorrect scatter
=#
for pokedf in [Age_p, Cas_p]
    pokedf.LogOut = log10.(pokedf.Out)
    pokedf.LogInterpoke = log10.(pokedf.PostInterpoke)
    transform!(pokedf, [:Streak, :Out, :LogOut] .=> zscore)
end
for streakdf in [Age_s, Cas_s]
    streakdf.LogDuration = log10.(streakdf.Trial_Duration)
    streakdf.ElapsedTime = log10.(streakdf.AfterLast_Duration .+ 0.1)
    streakdf.LogTravel = log10.(streakdf.Travel_to)
    streakdf.BinnedStreak = bin_axis(streakdf.Streak; unit_step = 4)
    streakdf.LongBinnedStreak = bin_axis(streakdf.Streak; unit_step = 20)
    streakdf.BlockBinnedStreak = bin_axis(streakdf.Streak; unit_step = 15)
    transform!(streakdf, :Streak .=> zscore)
end
open_html_table(FLPDevelopment.summarydf(Age_s,Age_p))
open_html_table(FLPDevelopment.summarydf(Cas_s,Cas_p))
########################### Example Session Plots ######################################
tt = filter(r -> r.MouseID == "RJ23" && 08<= r.Streak <= 52, Age_p)
plt = plot(legend = false, xlims = (-1,16), ylabel = "Trials",
    yaxis = false, xlabel = "Time (seconds)")
for r in eachrow(tt)
    FLPDevelopment.session_plot!(plt, r)
end
plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","AdultSession.pdf"))
##
tt = filter(r -> r.MouseID == "RJ02" && 08<= r.Streak <= 52, Age_p)
plt = plot(legend = false, xlims = (-1,16), ylabel = "Trials",
    yaxis = false, xlabel = "Time (seconds)")
for r in eachrow(tt)
    FLPDevelopment.session_plot!(plt, r)
end
plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","JuvenileSession.pdf"))
##################
comp_df = CSV.read("/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/Datasets/Stimulations/DRN_Opto_again/pokesDRN_Opto_again.csv", DataFrame)
comp90 = filter(r-> r.Protocol == "90/90", comp_df)
sessions = union(comp90.Session)
tt = filter(r-> r.Session == sessions[99] && 8<= r.Streak <= 52, comp90)
transform!(groupby(tt, [:Session, :Streak]), [:PokeIn, :PokeOut] => ((i,o) -> (In = i .- i[1], Out = o .- i[1])) => AsTable)
plt = plot(legend = false, xlims = (-1,16), ylabel = "Trials",
    yaxis = false, xlabel = "Time (seconds)")
for r in eachrow(tt)
    FLPDevelopment.session_plot!(plt, r)
end
plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","TrainedSession.pdf"))

########################### Over trials median Survival rate ######################################
Age_TrialMedianSurvival=plot();
att = combine(groupby(Age_s,:Age)) do dd
    mediansurvival_analysis(dd,:LogDuration,:LongBinnedStreak; plt = Age_TrialMedianSurvival)
end
Age_TrialMedianSurvival
xprop = ("Trial", (10,90), xyfont)
yprop = ("Median laeving time (log10 s)",(log10.([1,5,10,15,20]),[1,5,10,15,20]), xyfont)
plot!(xaxis = xprop, yaxis = yprop,
    margin = 1mm,
    left_margin = 10mm,
    size = (600,600), grid = true)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","2PartialMedianRatePerTrial.pdf"))
##
Cas_TrialMedianSurvival = plt = plot()
combine(groupby(Cas_s,:Virus)) do dd
    mediansurvival_analysis(dd,:LogDuration,:LongBinnedStreak; plt = Cas_TrialMedianSurvival)
end
Cas_TrialMedianSurvival
xaxis!(xlabel = "Trial", xlims = (10,90))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","PartialMedianRatePerTrial.pdf"))


###########################  Survival rate ######################################
Age_SurvivalAn, Age_SurvivalAn_df = function_analysis(Age_s,:LogDuration, survivalrate_algorythm; grouping = :Age, calc = :basic)
plot!(Age_SurvivalAn, xlabel = "Time (log10 s)", ylabel = "Survival rate", label = "")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","SurvivalRate.pdf"))
##
Cas_SurvivalAn, Cas_SurvivalAn_df = function_analysis(Cas_s,:LogDuration, survivalrate_algorythm; grouping = :Virus, calc = :basic)
plot!(Cas_SurvivalAn, xlabel = "Time (log10 s)", ylabel = "Survival rate", label = "")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","SurvivalRate.pdf"))


###########################  Model bootstrap ######################################
Age_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Basic = fit(MixedModel,Age_Basic_verb, Age_p, Bernoulli())
Age_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Age + LogOut_zscore * Age +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Full = fit(MixedModel,Age_Full_verb, Age_p, Bernoulli())
MixedModels.likelihoodratiotest(Age_Basic,Age_Full)
Age_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000AgeBootstrap.csv"),DataFrame)
# Age_BootDf = bootstrapdf(Age_p, Age_Full; n = 1000)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000AgeBootstrap.csv"),Age_BootDf)
leave_modelplt(Age_BootDf)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000ModelPlot.pdf"))
##
Cas_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +  (1+Streak_zscore+LogOut_zscore|MouseID));
Cas_Basic = fit(MixedModel,Cas_Basic_verb, Cas_p, Bernoulli())
Cas_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Virus + LogOut_zscore * Virus +  (1+Streak_zscore+LogOut_zscore|MouseID));
Cas_Full = fit(MixedModel,Cas_Full_verb, Cas_p, Bernoulli())
MixedModels.likelihoodratiotest(Cas_Basic,Cas_Full)
Cas_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","1000CasBootstrap.csv"), DataFrame)
# Cas_BootDf = bootstrapdf(Cas_p, Cas_Full, n = 1000)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","1000CasBootstrap.csv"),Cas_BootDf)
leave_modelplt(Cas_BootDf)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","1000ModelPlot.pdf"))
########################### Age Sex Model bootstrap ######################################
Age_Sex_verb = @formula(Leave ~ 1 + Streak_zscore*Sex + LogOut_zscore*Sex +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Sex = fit(MixedModel,Age_Sex_verb, Age_p, Bernoulli())
MixedModels.likelihoodratiotest(Age_Basic,Age_Sex)
# AgeSex_BootDf = bootstrapdf(Age_p, Age_Sex; n = 1000)
AgeSex_BootDf
leave_modelplt(AgeSex_BootDf; ylab = :names)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","1000SexModelPlot.pdf"))


########################### NumPokes scatter ######################################
Age_NP = Difference(Age_s, :Age, :Num_pokes; ylabel = "Median pokes per trial")
Age_NP.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","NP_test.pdf"))
##
Cas_NP = Difference(Cas_s, :Virus, :Num_pokes; ylabel = "Median pokes per trial")
Cas_NP.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","NP_test.pdf"))


########################## Incorrect scatter ######################################
Age_Inc = Difference(Age_s, :Age, :IncorrectLeave; ind_summary = incorrect_fraction, ylabel = "Fraction of premature leaving")
Age_Inc.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","Incorrect_test.pdf"))
##
Cas_Inc = Difference(Cas_s, :Virus, :IncorrectLeave; ind_summary = incorrect_fraction, ylabel = "Fraction of premature leaving")
Cas_Inc.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","Incorrect_test.pdf"))


########################### Travel scatter ######################################
Age_Trav = Difference(Age_s, :Age, :Travel_to; ind_summary = median, ylabel = "Median travel time (log10 s)")
Age_Trav.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","Travel_test.pdf"))
##
Cas_Trav = Difference(Cas_s, :Virus, :LogTravel; ind_summary = median, ylabel = "Median travel time (log10 s)")
Cas_Trav.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","Travel_test.pdf"))


########################### Interpoke scatter ######################################
Age_InterP = Difference(dropmissing(Age_p,:PostInterpoke), :Age, :PostInterpoke; ind_summary = median, ylabel = "Median inter-poke interval (log10 s)")
Age_InterP.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","Interpoke_test.pdf"))
##
Cas_InterP = Difference(dropmissing(Cas_p,:PostInterpoke), :Virus, :PostInterpoke; ind_summary = median, ylabel = "Median inter-poke interval (log10 s)")
Cas_InterP.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","Interpoke_test.pdf"))


########################### PokeDur scatter ######################################
Age_PDur = Difference(filter(r -> !r.Reward,dropmissing(Age_p,:PokeDur)), :Age, :PokeDur; ind_summary = median, ylabel = "Median poke duration (s)")
Age_PDur.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","pokeDur_test.pdf"))
##
Cas_PDur = Difference(filter(r -> !r.Reward,dropmissing(Cas_p,:PokeDur)), :Virus, :PokeDur; ind_summary = median, ylabel = "Median poke duration (s)")
Cas_PDur.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","pokeDur_test.pdf"))


########################### Duty cycle ######################################
Age_df = combine(groupby(Age_p, [:MouseID,:Streak, :Age]), :PokeDur => sum => :PokingTime, [:PokeIn, :PokeOut] => ((i,o) -> o[end] - i[1]) => :TrialDuration)
Age_df.RelativePokingTime = Age_df.PokingTime./Age_df.TrialDuration
Age_Poking = Difference(Age_df, :Age, :RelativePokingTime; ind_summary = median, ylabel = "Median fraction of time spent poking during trial")
Age_Poking.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","DutyCycle_test.pdf"))
###
Cas_df = combine(groupby(Cas_p, [:MouseID,:Streak, :Virus]), :PokeDur => sum => :PokingTime, [:PokeIn, :PokeOut] => ((i,o) -> o[end] - i[1]) => :TrialDuration)
Cas_df.RelativePokingTime = Cas_df.PokingTime ./ Cas_df.TrialDuration
Cas_Poking = Difference(Cas_df, :Virus, :RelativePokingTime; ind_summary = median, ylabel = "Median fraction of time spent poking during trial")
Cas_Poking.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","DutyCycle_test.pdf"))



#######################
Age_Full = fit(MixedModel,Age_Full_verb, Age_p, Bernoulli())
Cas_Full = fit(MixedModel,Cas_Full_verb, Cas_p, Bernoulli())

Age_p.Model = predict(Age_Full)
Cas_p.Model = predict(Cas_Full)
Age_mcheck = Age_p[:,[:MouseID, :Age, :Streak_zscore, :LogOut_zscore, :Leave, :Model]]
Cas_mcheck = Cas_p[:,[:MouseID, :Virus, :Streak_zscore, :LogOut_zscore, :Leave, :Model]]

simulate!(Age_Full)
Age_Full

#######################################################EXCLUDED PLOTS #########################################
########################### Median Survival rate ######################################
Age_MedianSurvival = mediansurvival_analysis(Age_s,:LogDuration, :Age)
xaxis!(xlims = (-0.25,2.25), xlabel = "Group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","MedianRate.pdf"))
####
Cas_MedianSurvival = mediansurvival_analysis(Cas_s,:LogDuration, :Virus)
xaxis!(xlims = (-0.25,2.25), xlabel = "Group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","MedianRate.pdf"))
########################### Age Survival rate per trial ######################################
agedf = filter(r -> 10 <= r.Streak <= 70, Age_s)
agedf.BlockBinnedStreak = bin_axis(agedf.Streak; unit_step = 30)
Age_SurvivalAn = function_analysis(agedf,:LogDuration, survivalrate_algorythm; grouping = :BlockBinnedStreak,
    color = [:burlywood4], linestyle = [:solid :dash])
plot!(ylabel = "SurvivalRate", xlabel = "Time (log10 s)", legendtitle = "Trial group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","Filt10-70_TrialSurvivalRate.pdf"))
agedf = filter(r -> r.Streak <= 71, Age_s)
agedf.BlockBinnedStreak = bin_axis(agedf.Streak; unit_step = 35)
Age_SurvivalAn = function_analysis(agedf,:LogDuration, survivalrate_algorythm; grouping = :BlockBinnedStreak,
    color = [:burlywood4], linestyle = [:solid :dash])
plot!(ylabel = "SurvivalRate", xlabel = "Time (log10 s)", legendtitle = "Trial group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","TrialSurvivalRate.pdf"))
########################### Cas Survival rate per trial ######################################
casdf = filter(r -> 10 <= r.Streak <= 70, Cas_s)
casdf.BlockBinnedStreak = bin_axis(casdf.Streak; unit_step = 30)
Cas_SurvivalAn = function_analysis(casdf,:LogDuration, survivalrate_algorythm; grouping = :BlockBinnedStreak,
    color = [:burlywood4], linestyle = [:solid :dash])
plot!(ylabel = "SurvivalRate", xlabel = "Time (log10 s)", legendtitle = "Trial group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","Filt10-70_TrialSurvivalRate.pdf"))
casdf = filter(r -> r.Streak <= 71, Cas_s)
casdf.BlockBinnedStreak = bin_axis(casdf.Streak; unit_step = 35)
Cas_SurvivalAn = function_analysis(casdf,:LogDuration, survivalrate_algorythm; grouping = :BlockBinnedStreak,
    color = [:burlywood4], linestyle = [:solid :dash])
plot!(ylabel = "SurvivalRate", xlabel = "Time (log10 s)", legendtitle = "Trial group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","TrialSurvivalRate.pdf"))
########################### Age Leaving time over trial ######################################
Age_OverTrial_df = summary_xy(Age_s,:BinnedStreak,:LogDuration; group = :Age)
sort!(Age_OverTrial_df,[:Age,:BinnedStreak])
@df Age_OverTrial_df plot(string.(:BinnedStreak),:Mean, group = :Age, ribbon = :Sem,
    linecolor = :auto, xlabel = "Trial", ylabel = "Leaving time (log10 s)",size=(650,600))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","LeavingOverTrial.pdf"))
########################### Virus Leaving time over trial ######################################
Cas_OverTrial_df = summary_xy(Cas_s,:BinnedStreak,:LogDuration; group = :Virus)
sort!(Cas_OverTrial_df,[:Virus,:BinnedStreak])
@df Cas_OverTrial_df plot(string.(:BinnedStreak),:Mean, group = :Virus, ribbon = :Sem,
    linecolor = :auto, xlabel = "Trial", ylabel = "Leaving time (log10 s)",size=(650,600))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","LeavingOverTrial.pdf"))
