## Modules call
using Revise, FLPDevelopment, BrowseTables
## Plot settings
# Font size and type
xyfont = font(18, "Bookman Light")
legfont = font(14, "Bookman Light")

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
    guidefont = xyfont,
    legendfont = legfont)
## Load and adjust data
include("Young_to_run2.jl")
include("Caspase_to_run.jl")
for df in (Age_p, Age_b, Age_s, Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Protocol == "90/90" &&
    r.MouseID != "CD09" && # biting, see B1_CD09_2020-07-13 minute30
    r.MouseID != "RJ58" && # blind
    r.MouseID != "RJ67" && # biting, see B3_RJ67_2020-09-28 minute 7:33
    !(r.MouseID in first_females_group) &&
    r.ProtocolSession == 1
    ,df)
end
for df in (Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Gen == "Rbp4-cre", df) # exclude wild type animals
end
#adjustments for pokes DataFrames
for pokedf in [Age_p, Cas_p]
    # transform time in log 10 scale
    pokedf.LogOut = log10.(pokedf.Out)
    pokedf.LogInterpoke = log10.(pokedf.PostInterpoke)
    # zscore values for logistic regression
    transform!(pokedf, [:Streak, :Out, :LogOut] .=> zscore)
end
#adjustments for trials DataFrames
for streakdf in [Age_s, Cas_s]
    # transform time in log 10 scale
    streakdf.LogDuration = log10.(streakdf.Trial_Duration)
    streakdf.ElapsedTime = log10.(streakdf.AfterLast_Duration .+ 0.1)
    streakdf.LogTravel = log10.(streakdf.Travel_to)
    # bin trials
    # streakdf.BinnedStreak = bin_axis(streakdf.Streak; unit_step = 4)
    # streakdf.BlockBinnedStreak = bin_axis(streakdf.Streak; unit_step = 15)
    streakdf.LongBinnedStreak = bin_axis(streakdf.Streak; unit_step = 20)
    transform!(streakdf, :Streak .=> zscore)
end
open_html_table(FLPDevelopment.summarydf(Age_s,Age_p))
open_html_table(FLPDevelopment.summarydf(Cas_s,Cas_p))
########################### Example Session Plots ######################################
## Adult example session
lowertrial = 20
uppertrial = 45
tt = filter(r -> r.MouseID == "RJ23" && lowertrial <= r.Streak <= uppertrial, Age_p)
session_plot(tt)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","AdultSession.pdf"))
## Juvenile example session
tt = filter(r -> r.MouseID == "RJ02" && lowertrial <= r.Streak <= uppertrial, Age_p)
session_plot(tt)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","JuvenileSession.pdf"))
## Well trained example session
comp_df = CSV.read("/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/Datasets/Stimulations/DRN_Opto_again/pokesDRN_Opto_again.csv", DataFrame)
comp90 = filter(r-> r.Protocol == "90/90", comp_df)
sessions = union(comp90.Session)
tt = filter(r-> r.Session == sessions[99] && lowertrial <= r.Streak <= uppertrial, comp90)
transform!(groupby(tt, [:Session, :Streak]), [:PokeIn, :PokeOut] => ((i,o) -> (In = i .- i[1], Out = o .- i[1])) => AsTable)
session_plot(tt)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","TrainedSession.pdf"))
########################### Poking rate ######################################
# Age experiment
# calculate the total amount of time spent in poking and the trial duration
Age_df = combine(groupby(Age_p, [:MouseID,:Streak, :Age]), :PokeDur => sum => :PokingTime,
    [:PokeIn, :PokeOut] => ((i,o) -> o[end] - i[1]) => :TrialDuration)
# calculate the fraction of time spent poking per trial: Poking rate
Age_df.PokingRate = Age_df.PokingTime./Age_df.TrialDuration
# Test the difference in the poking rate for age group
Age_Poking = Difference(Age_df, :Age, :PokingRate; ind_summary = median, ylabel = "Poking rate", ylims = (0,1))
Age_Poking.plt
Age_Poking.groupdf
Age_Poking.test
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","DutyCycle_test.pdf"))
## Same for ablation experiment
Cas_df = combine(groupby(Cas_p, [:MouseID,:Streak, :Virus]), :PokeDur => sum => :PokingTime, [:PokeIn, :PokeOut] => ((i,o) -> o[end] - i[1]) => :TrialDuration)
Cas_df.PokingRate = Cas_df.PokingTime ./ Cas_df.TrialDuration
Cas_Poking = Difference(Cas_df, :Virus, :PokingRate; ind_summary = median, ylabel = "Poking rate", ylims = (0,1))
Cas_Poking.plt
Cas_Poking.groupdf
Cas_Poking.test
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","DutyCycle_test.pdf"))
###########################  Leving rate ######################################
Age_LevingAn, Age_LevingAn_df = function_analysis(Age_s,:LogDuration, cumulative_algorythm; grouping = :Age, calc = :bootstrapping)
xprop = ("Poke Time (log10 s)", xyfont,(log10.([0.1,1,10,100,1000]),["0.1","1","10","100","1000"]))
yprop = ("Probablity of leaving", xyfont)
plot!(Age_LevingAn, xaxis = xprop, yaxis = yprop, legend = false)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","LevingRate.pdf"))
##
Cas_SurvivalAn, Cas_SurvivalAn_df = function_analysis(Cas_s,:LogDuration, cumulative_algorythm; grouping = :Virus, calc = :bootstrapping)
plot!(Cas_SurvivalAn, xaxis = xprop, yaxis = yprop, legend = false)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","LevingRate.pdf"))
###########################  Model bootstrap ######################################
Age_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Basic = fit(MixedModel,Age_Basic_verb, Age_p, Bernoulli())
Age_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Age + LogOut_zscore * Age +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Full = fit(MixedModel,Age_Full_verb, Age_p, Bernoulli())
MixedModels.likelihoodratiotest(Age_Basic,Age_Full)
Age_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000AgeBootstrap.csv"),DataFrame)
# Age_BootDf = bootstrapdf(Age_p, Age_Full; n = 1000)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000AgeBootstrap.csv"),Age_BootDf)
Age_btdf = Age_BootDf[[1,4,2,3,5,6],:]
Age_btdf.err = [tuple(parse.(Float64,split(err[2:end-1],", "))...) for err in Age_btdf.err]
Age_btdf.variable
yprop = ("",font(10, "Bookman Light"),(collect(1:nrow(Age_btdf)),
    [L"Intercept",L"PokeTime",L"Trial",L"Juveniles",
    L"Trial \&",
    L"PokeTime \&"]))
    xprop = ("Coefficient estimate", xyfont, (-1.33,1.3))
    Titolo = L"+ \Leftarrow Behavioural Control \Rightarrow -"
    @df Age_btdf scatter(:coef ,1:nrow(Age_btdf), xerror = :err,legend = false,
    xaxis = xprop, yaxis = yprop, markercolor = :gray75, title = Titolo)
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000AgeModelPlot.pdf"))
# leave_modelplt(Age_BootDf)
##
Cas_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +  (1+Streak_zscore+LogOut_zscore|MouseID));
Cas_Basic = fit(MixedModel,Cas_Basic_verb, Cas_p, Bernoulli())
Cas_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Virus + LogOut_zscore * Virus +  (1+Streak_zscore+LogOut_zscore|MouseID));
Cas_Full = fit(MixedModel,Cas_Full_verb, Cas_p, Bernoulli())
MixedModels.likelihoodratiotest(Cas_Basic,Cas_Full)
# Cas_BootDf = bootstrapdf(Cas_p, Cas_Full, n = 1000)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","1000CasBootstrap.csv"),Cas_BootDf)
Cas_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","1000CasBootstrap.csv"), DataFrame)
Cas_btdf = Cas_BootDf[[1,4,2,3,5,6],:]
Cas_btdf.err = [tuple(parse.(Float64,split(err[2:end-1],", "))...) for err in Cas_btdf.err]
Cas_btdf.variable
yprop = ("",font(10, "Bookman Light"),(collect(1:nrow(Cas_btdf)),
    [L"Intercept",L"PokeTime",L"Trial",L"Caspase",
    L"Trial \&",
    L"PokeTime\&"]))
    xprop = ("Coefficient estimate", xyfont, (-1.4,1.35))
    @df Cas_btdf scatter(:coef ,1:nrow(Cas_btdf), xerror = :err,legend = false,
    xaxis = xprop, yaxis = yprop, markercolor = :gray75, title = Titolo)
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","1000ModelPlot.pdf"))
# leave_modelplt(Cas_BootDf)


########################## Incorrect scatter ######################################
Age_Inc = Difference(Age_s, :Age, :IncorrectLeave; ind_summary = incorrect_fraction, ylabel = "Fraction of premature leaving", ylims = (0,1))
Age_Inc.plt
Age_Inc.groupdf
Age_Inc.test
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","Incorrect_test.pdf"))
##
Cas_Inc = Difference(Cas_s, :Virus, :IncorrectLeave; ind_summary = incorrect_fraction, ylabel = "Fraction of premature leaving", ylims = (0,1))
Cas_Inc.plt
Cas_Inc.groupdf
Cas_Inc.test
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","Incorrect_test.pdf"))
########################### Age Sex Model bootstrap ######################################
###########################  Leving rate ######################################
Sex_LevingAn, Sex_LevingAn_df = function_analysis(Age_s,:LogDuration, cumulative_algorythm;
    grouping = :Sex, calc = :bootstrapping, color = [:violetred3 :sandybrown])
    xprop = ("Poke Time (log10 s)", xyfont,(log10.([0.1,1,10,100,1000]),["0.1","1","10","100","1000"]))
    yprop = ("Probablity of leaving", xyfont)
    plot!(Sex_LevingAn, xaxis = xprop, yaxis = yprop, legend = (0.8, 0.25))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","SexLevingRate.pdf"))

########################### Age Sex Model bootstrap #####################################
Age_Sex_verb1 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + Age +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Sex_verb2 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + Age*Sex +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Sex1 = fit(MixedModel,Age_Sex_verb1, Age_p, Bernoulli())
Age_Sex2 = fit(MixedModel,Age_Sex_verb2, Age_p, Bernoulli())
MixedModels.likelihoodratiotest(Age_Sex1,Age_Sex2)
# Sex_BootDf2 = bootstrapdf(Age_p, Age_Sex2; n = 1000)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","1000SexBootstrap.csv"),Sex_BootDf2)
Sex_BootDf2 = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","1000SexBootstrap.csv"),DataFrame)
Sex_btdf = Sex_BootDf2[[1,3,2,4,5,6],:]
Sex_btdf.err = [tuple(parse.(Float64,split(err[2:end-1],", "))...) for err in Sex_btdf.err]
yprop = ("",font(10, "Bookman Light"),(collect(1:nrow(Sex_btdf)),
    [L"Intercept",L"PokeTime",L"Trial",L"Juveniles",
    L"Male",
    L"Juveniles \&"]))
    xprop = ("Coefficient estimate", xyfont, (-1.33,1.3))
    @df Sex_btdf scatter(:coef ,1:nrow(Sex_btdf), xerror = :err,legend = false,
    xaxis = xprop, yaxis = yprop, markercolor = :gray75, title = Titolo)
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","1000FullSexModelPlot.pdf"))
########################### Over trials median Survival rate ######################################
Age_TrialMedianSurvival=plot();
att = combine(groupby(Age_s,:Age)) do dd
    mediansurvival_analysis(dd,:LogDuration,:LongBinnedStreak; plt = Age_TrialMedianSurvival)
end
Age_TrialMedianSurvival
xprop = ("Trial", (10,90), xyfont)
yprop = ("Median laeving time (log10 s)", xyfont,(log(1),log10(30)),(log10.([1,5,10,15,20, 30]),string.([1,5,10,15,20,30])))
plot!(xaxis = xprop, yaxis = yprop, grid = true, size = (600,600))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","PartialMedianRatePerTrial.pdf"))
##
Cas_TrialMedianSurvival = plt = plot()
combine(groupby(Cas_s,:Virus)) do dd
    mediansurvival_analysis(dd,:LogDuration,:LongBinnedStreak; plt = Cas_TrialMedianSurvival)
end
Cas_TrialMedianSurvival
plot!(xaxis = xprop, yaxis = yprop, grid = true, size = (600,600))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","PartialMedianRatePerTrial.pdf"))
