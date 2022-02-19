using Revise, FLPDevelopment, BrowseTables
# Plot settings
# Font size and type
xyfont = font(18, "Bookman Light")
legfont = font(14, "Bookman Light")

# Figure/Plot image size
pixel_factor = 600
ps_x = pixel_factor * 1
ps_y = pixel_factor * 1
fig_size = (ps_x, ps_y)

theme(:default)
# gr(size = fig_size,
#     titlefont = font(14, "Bookman Light"),
#     guidefont = font(14, "Bookman Light"),
#     tick_orientation = :out,
#     grid = false,
#     markerstrokecolor = :black,
#     markersize = 8,
#     thickness_scaling = 1.5,
#     legendfont = font(10, "Bookman Light"))
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
    pokedf.LogDuration = log10.(pokedf.PokeDur)
    # zscore values for logistic regression
    transform!(pokedf, [:Streak, :Out, :LogOut] .=> zscore)
    transform!(groupby(pokedf,[:MouseID, :Streak]), :PokeDur => cumsum => :CumDuration)
    pokedf.Occupancy = pokedf.CumDuration ./ pokedf.Out
    pokedf.Occupancy_zscore = zscore(pokedf.Occupancy)
end
#adjustments for trials DataFrames
for streakdf in [Age_s, Cas_s]
    # transform time in log 10 scale
    streakdf.LogDuration = log10.(streakdf.Trial_Duration)
    # streakdf.ElapsedTime = log10.(streakdf.AfterLast_Duration .+ 0.1)
    streakdf.LogTravel = log10.(streakdf.Travel_to)
    # bin trials
    # streakdf.BinnedStreak = bin_axis(streakdf.Streak; unit_step = 4)
    # streakdf.BlockBinnedStreak = bin_axis(streakdf.Streak; unit_step = 15)
    streakdf.LongBinnedStreak = bin_axis(streakdf.Streak; unit_step = 20)
    transform!(streakdf, :Streak .=> zscore)
end
open_html_table(FLPDevelopment.summarydf(Age_s,Age_p))
open_html_table(FLPDevelopment.summarydf(Cas_s,Cas_p))
###

for df in [Age_p, Age_s, Cas_p, Cas_s]
    x = :Age in propertynames(df) ? :Age : :Virus
    caso = :Age in propertynames(df) ? "Age" : "Virus"
    y = :Correct in propertynames(df) ? [:MouseID, x,:Sex,:Day,:Side, :Streak,:Poke,:Correct,:Reward,:In,:Out] :
        [:MouseID,x,:Sex,:Day,:Side, :Streak,:IncorrectLeave, :Num_pokes, :Num_Rewards, :Trial_Duration]
    tipo = :Correct in propertynames(df) ? "Pokes" : "Trials"
    prov = select(df,y)
    rename!(prov, :Streak => :Trial)
    :Num_pokes in propertynames(prov) && rename!(prov, :Num_pokes => :Num_Pokes)
    fname = joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","FullData" * tipo * caso * ".csv")
    CSV.write(fname,prov)
end
########################### Fig2 ######################################
########################### B ###########################
########################### Example Session Plots ######################################
## Adult example session
lowertrial = 20
uppertrial = 45
# ex_ad = filter(r -> r.MouseID == "RJ23" && lowertrial <= r.Streak <= uppertrial, Age_p)
ex_ad = filter(r -> r.MouseID == "RJ23", Age_p)
t1 = filter(r -> r.Leave, ex_ad)
sort!(t1,[:Out])
order = Dict(t1.Streak .=> 1:nrow(t1))
ex_ad.order = [get(order,x,0) for x in ex_ad.Streak]
sort!(ex_ad, :order)
filter!(r -> r.order != 0, ex_ad)
plt = session_plot(ex_ad; ordered = true)
xaxis!(plt, xlims = :auto)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "B", "AdultSession.pdf"))
select!(ex_ad,[:Day,:MouseID,:Age,:order,:In,:Out,:Reward,:Correct])
## Juvenile example session
# ex_juv = filter(r -> r.MouseID == "RJ02" && lowertrial <= r.Streak <= uppertrial, Age_p)
ex_juv = filter(r -> r.MouseID == "RJ02", Age_p)
t1 = filter(r -> r.Leave, ex_juv)
sort!(t1,[:Out])
order = Dict(t1.Streak .=> 1:nrow(t1))
ex_juv.order = [get(order,x,0) for x in ex_juv.Streak]
sort!(ex_juv, :order)
filter!(r -> r.order != 0, ex_juv)
plt = session_plot(ex_juv; ordered = true)
xaxis!(plt, xlims = :auto)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "B", "JuvenileSession.pdf"))
select!(ex_juv,[:Day,:MouseID,:Age,:order,:In,:Out,:Reward,:Correct])
## Well trained example session
# comp_df = CSV.read("/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/Datasets/Stimulations/DRN_Opto_again/pokesDRN_Opto_again.csv", DataFrame)
# comp_df = DataFrame(CSV.File(joinpath("/Volumes/GoogleDrive/My Drive/Flipping/Datasets/","Stimulations/DRN_Opto_again/pokesDRN_Opto_again.csv")))
mother1 = replace(path,basename(path)=>"")[1:end-1]
mother2 = replace(mother1,basename(mother1)=>"")
comp_df = CSV.read(joinpath(mother2,"Stimulations","DRN_Opto_again","pokesDRN_Opto_again.csv"), DataFrame)
comp90 = filter(r-> r.Protocol == "90/90", comp_df)
sessions = union(comp90.Session)
# ex_train = filter(r-> r.Session == sessions[99] && lowertrial <= r.Streak <= uppertrial, comp90)
ex_train = filter(r-> r.Session == sessions[99], comp90)
transform!(groupby(ex_train, [:Session, :Streak]), [:PokeIn, :PokeOut] => ((i,o) -> (In = i .- i[1], Out = o .- i[1])) => AsTable)
filter!(r-> r.Leave != "NA", ex_train)
t1 = filter(r -> parse(Bool,r.Leave), ex_train)
sort!(t1,[:Out])
order = Dict(t1.Streak .=> 1:nrow(t1))
ex_train.order = [get(order,x,0) for x in ex_train.Streak]
sort!(ex_train, :order)
filter!(r -> r.order != 0, ex_train)
session_plot(ex_train; ordered = true)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "B", "TrainedSession.pdf"))
ex_train[!,:Age] .= "Trained"
select!(ex_train,[:Day,:MouseID,:Age,:order,:In,:Out,:Reward,:Correct])
tt = vcat(ex_juv,ex_ad,ex_train)
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "B","PokeData.csv"), ex_juv)
########################### C ###########################
## Density leaving times
@df Age_s density(:LogDuration, group = :Age,
    xlabel = "Leaving time (s)", ylabel = "Kernell density estimate", xticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "C", "DensityTrialDuration.pdf"))
Ad = kde(Age_s[Age_s.Age .== "Adults",:LogDuration])
Juv = kde(Age_s[Age_s.Age .== "Juveniles",:LogDuration])
prov = vcat(DataFrame(LeavingTime = collect(Ad.x), KernelDensity = Ad.density, Age = "Adults"),
    DataFrame(LeavingTime = collect(Juv.x), KernelDensity = Juv.density, Age = "Juveniles"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "C","TrialDurationData.csv"), prov)
# ########################## D ###########################
## Age Survival Plot
Age_LevingAn, Age_LevingAn_df = function_analysis(Age_s,:LogDuration, cumulative_algorythm; grouping = :Age, calc = :bootstrapping)
xprop = ("Poke Time(seconds)", xyfont,(log10.([0.1,1,10,100,1000]),["0.1","1","10","100","1000"]))
yprop = ("Probablity of leaving", xyfont)
plot!(Age_LevingAn, xaxis = xprop, yaxis = yprop, legend = false)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2","E","LevingRate.pdf"))
rename!(Age_LevingAn_df, :LogDuration => :PokeTime)
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "D","LevingRate.csv"), Age_LevingAn_df)
########################### E ###########################
##
Age_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Basic = fit(MixedModel,Age_Basic_verb, Age_p, Bernoulli())
Age_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Age + LogOut_zscore * Age +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Full = fit(MixedModel,Age_Full_verb, Age_p, Bernoulli())
MixedModels.likelihoodratiotest(Age_Basic,Age_Full)
Age_noTrial_verb = @formula(Leave ~ 1 + LogOut_zscore * Age +  (1+LogOut_zscore|MouseID));
Age_noTrial =  fit(MixedModel,Age_noTrial_verb, Age_p, Bernoulli())
MixedModels.likelihoodratiotest(Age_noTrial,Age_Full)
Age_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000AgeBootstrap.csv"),DataFrame)
# Age_BootDf = bootstrapdf(Age_p, Age_Full; n = 1000)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000AgeBootstrap.csv"),Age_BootDf)
Age_btdf = Age_BootDf[[1,4,2,3,5,6],:]
Age_btdf.err = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in Age_btdf.err]
# Age_btdf.variable
yprop = ("",font(10, "Bookman Light"),(collect(1:nrow(Age_btdf)),
    [L"Intercept",L"PokeTime",L"Trial",L"Juveniles",
    L"Trial \&",
    L"PokeTime \&"]))
    xprop = ("Coefficient estimate", xyfont, (-1.33,1.3))
    Titolo = L"+ \Leftarrow Behavioural Control \Rightarrow -"
    @df Age_btdf scatter(:coef ,1:nrow(Age_btdf), xerror = :err,legend = false,
    xaxis = xprop, yaxis = yprop, markercolor = :gray75)
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2","E","LevingRateBootstrap.pdf"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "E","LevingRateBootstrap.csv"),
    select(Age_btdf, Not([:type,:group])))
########################## F1 ###########################
Age_Npokes = Difference(Age_s, :Age, :Num_pokes, ylabel = "number of pokes per trial", ylims = (0,5.5))
Age_Npokes.plt
Age_Npokes.test
Age_Npokes.groupdf
Age_Npokes.groupdf[!,:Measure] .= "Number of pokes"
rename!(Age_Npokes.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "F","NumPokes_test.pdf"))
########################## F2 ###########################
Age_Inc = Difference(Age_s, :Age, :IncorrectLeave; ind_summary = incorrect_fraction, ylabel = "Fraction of premature leaving", ylims = (0,1))
Age_Inc.plt
Age_Inc.test
Age_Inc.groupdf
Age_Inc.groupdf[!,:Measure] .= "Fraction of premature leaving"
rename!(Age_Inc.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "F","Incorrect_test.pdf"))
prov = vcat(Age_Npokes.groupdf,Age_Inc.groupdf)
select!(prov, Not(:color))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "F","Tests.csv"), prov)
########################### G ###########################
##Port occupancy by pokes
# df0 = FLPDevelopment.summarizexy(Age_p,:LogOut,:Occupancy, group = :Age, bin = true, digits = 1)
# @df filter(r -> !isnan(r.Occupancy_sem), df0) plot(:BinnedLogOut,:Occupancy_mean, ribbon = :Occupancy_sem, group = :Age,
#     xlabel = "Elapsed time (Log10 seconds)", xticks = -1:3, ylabel = "Fraction of time poking", legend = :topright)
AgeOcc = FLPDevelopment.summarizexy(Age_p,:LogOut,:Occupancy, group = :Age, bin = true, digits = 1, calc = :bootstrapping)
rename!(AgeOcc,[:Central => :Median, :ERR => :CI])
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "G","AgeOccupancy.csv"), select(AgeOcc, Not([:low,:up])))
@df AgeOcc plot(:BinnedLogOut,:Median, ribbon = (:low, :up), group = :Age,
    xlabel = "Elapsed time (Log10 seconds)", xticks = -1:3, ylabel = "Fraction of time poking", legend = :topright)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2","G","AgeOccupancy.pdf"))
########################### H ###########################
## Port occupancy model
occupancy_df = Age_p
occupancy_df.Occupancy_zscore = zscore(occupancy_df.Occupancy)
Occupancy_v1 = @formula(Occupancy_zscore ~ 1 + LogOut_zscore + (1+LogOut_zscore|MouseID));
Occupancy_v2 = @formula(Occupancy_zscore ~ 1 + LogOut_zscore * Age + (1+LogOut_zscore|MouseID));
Occupancy_m1 = fit(MixedModel,Occupancy_v1, occupancy_df)
Occupancy_m2 = fit(MixedModel,Occupancy_v2, occupancy_df)
MixedModels.likelihoodratiotest(Occupancy_m1,Occupancy_m2)
# OccupancyAge_BootDf = bootstrapdf(occupancy_df, Occupancy_m2; n = 1000)
# OccupancyAge_BootDf.names = ["Intercept", "PokeTime","Group:Juveniles","PokeTime&Group:Juveniles"]
# OccupancyAge_BootDf.Y = 1:nrow(OccupancyAge_BootDf)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","ZMfeedback","Occupancy","1000AgeOccupancyBootstrap.csv"),OccupancyAge_BootDf)
OccupancyAge_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","ZMfeedback","Occupancy","1000AgeOccupancyBootstrap.csv"), DataFrame)
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2", "H","AgeOccupancyBootstrap.csv"),
    select(OccupancyAge_BootDf, Not([:type,:group,:Y])))
OccupancyAge_BootDf[!,:err] = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in OccupancyAge_BootDf.err]
xprop = (L"Coefficient estimate", xyfont, (-1.1,0.2))
@df OccupancyAge_BootDf scatter(:coef ,1:nrow(OccupancyAge_BootDf), xerror = :err,legend = false,
    markercolor = :gray75, size = (600,600),
    xaxis = xprop,
    yticks = (collect(1:4), [L"Intercept",L"PokeTime",L"Group:Juveniles",L"PokeTime \& Group:Juveniles"]))
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
    annotate!([(-1,1,"n.s", 10), (-1,2,"*", 10), (-1,3,"*", 10), (-1,4,"n.s", 10)])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig2","H","AgeOccupancyBootstrap.pdf"))
########################### Fig3 ###########################
########################### E ###########################
@df Cas_s density(:LogDuration, group = :Virus,
    xlabel = "Leaving time (s)", ylabel = "Kernell density estimate", xticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "E", "DensityTrialDuration.pdf"))
Cs = kde(Cas_s[Cas_s.Virus .== "Caspase",:LogDuration])
Td = kde(Cas_s[Cas_s.Virus .== "tdTomato",:LogDuration])
prov = vcat(DataFrame(LeavingTime = collect(Cs.x), KernelDensity = Cs.density, Virus = "Caspase"),
    DataFrame(LeavingTime = collect(Td.x), KernelDensity = Td.density, Virus = "tdTomato"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "E","TrialDurationData.csv"), prov)
########################## F1 ###########################
Cas_Npokes = Difference(Cas_s, :Virus, :Num_pokes, ylabel = "number of pokes per trial", ylims = (0,5.5))
Cas_Npokes.plt
Cas_Npokes.test
Cas_Npokes.groupdf
Cas_Npokes.groupdf[!,:Measure] .= "Number of pokes"
rename!(Cas_Npokes.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "F","NumPokes_test.pdf"))
########################## F2 ###########################
Cas_Inc = Difference(Cas_s, :Virus, :IncorrectLeave; ind_summary = incorrect_fraction, ylabel = "Fraction of premature leaving", ylims = (0,1))
Cas_Inc.plt
Cas_Inc.test
Cas_Inc.groupdf
Cas_Inc.groupdf[!,:Measure] .= "Fraction of premature leaving"
rename!(Cas_Inc.groupdf,[:Central => :Median, :ERR => :CI])
prov = vcat(Cas_Npokes.groupdf,Cas_Inc.groupdf)
select!(prov, Not(:color))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "F","Incorrect_test.pdf"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "F","Tests.csv"), prov)
########################### G ###########################
##
Cas_SurvivalAn, Cas_SurvivalAn_df = function_analysis(Cas_s,:LogDuration, cumulative_algorythm; grouping = :Virus, calc = :bootstrapping)
xprop = ("Poke Time(seconds)", xyfont,(log10.([0.1,1,10,100,1000]),["0.1","1","10","100","1000"]))
yprop = ("Probablity of leaving", xyfont)
plot!(Cas_SurvivalAn, xaxis = xprop, yaxis = yprop, legend = false)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3","G","VirusLeavingRate.pdf"))
rename!(Cas_SurvivalAn_df, [:Central => :Median, :ERR => :CI, :LogDuration => :PokeTime])
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "G","VirusLevingRate.csv"), Cas_SurvivalAn_df)
########################### H ###########################
##
Cas_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +  (1+Streak_zscore+LogOut_zscore|MouseID));
Cas_Basic = fit(MixedModel,Cas_Basic_verb, Cas_p, Bernoulli())
Cas_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Virus + LogOut_zscore * Virus +  (1+Streak_zscore+LogOut_zscore|MouseID));
Cas_Full = fit(MixedModel,Cas_Full_verb, Cas_p, Bernoulli())
MixedModels.likelihoodratiotest(Cas_Basic,Cas_Full)
# Cas_BootDf = bootstrapdf(Cas_p, Cas_Full, n = 1000)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","1000CasBootstrap.csv"),Cas_BootDf)
Cas_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","1000CasBootstrap.csv"), DataFrame)
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "H","VirusLevingRatBootstrap.csv"),
    select(Cas_BootDf, Not([:type,:group])))
Cas_btdf = Cas_BootDf[[1,4,2,3,5,6],:]
Cas_btdf.err = [tuple(parse.(Float64,split(err[2:end-1],", "))...) for err in Cas_btdf.err]
Cas_btdf.variable
yprop = ("",font(10, "Bookman Light"),(collect(1:nrow(Cas_btdf)),
    [L"Intercept",L"PokeTime",L"Trial",L"Caspase",
    L"Trial \&",
    L"PokeTime\&"]))
    xprop = ("Coefficient estimate", xyfont, (-1.4,1.35))
    @df Cas_btdf scatter(:coef ,1:nrow(Cas_btdf), xerror = :err,legend = false,
    xaxis = xprop, yaxis = yprop, markercolor = :gray75)
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "H","VirusLevingRatBootstrap.pdf"))
########################### I ###########################
##Port occupancy by pokes
# df = FLPDevelopment.summarizexy(Cas_p,:LogOut,:Occupancy, group = :Virus, bin = true, digits = 1)
# filter!(r -> !isnan(r.Occupancy_sem), df)
# @df df plot(:BinnedLogOut,:Occupancy_mean, ribbon = :Occupancy_sem, group = :Virus,
#     xlabel = "Elapsed time (Log10 seconds)", xticks = -1:3, ylabel = "Fraction of time poking", legend = :topright)
CasOcc = FLPDevelopment.summarizexy(Cas_p,:LogOut,:Occupancy, group = :Virus, bin = true, digits = 1, calc = :bootstrapping)
rename!(CasOcc,[:Central => :Median, :ERR => :CI])
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "I","VirusOccupancy.csv"), select(CasOcc, Not([:low,:up])))
@df CasOcc plot(:BinnedLogOut,:Median, ribbon = (:low, :up), group = :Virus,
    xlabel = "Elapsed time (Log10 seconds)", xticks = -1:3, ylabel = "Fraction of time poking", legend = :topright)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "I","VirusOccupancy.pdf"))
########################### J ###########################
## Port occupancy model
occupancy_df = Cas_p
Occupancy_v1 = @formula(Occupancy_zscore ~ 1 + LogOut_zscore + (1+LogOut_zscore|MouseID));
Occupancy_v2 = @formula(Occupancy_zscore ~ 1 + LogOut_zscore * Virus + (1+LogOut_zscore|MouseID));
Occupancy_m1 = fit(MixedModel,Occupancy_v1, occupancy_df)
Occupancy_m2 = fit(MixedModel,Occupancy_v2, occupancy_df)
MixedModels.likelihoodratiotest(Occupancy_m1,Occupancy_m2)
# OccupancyCas_BootDf = bootstrapdf(occupancy_df, Occupancy_m2; n = 1000)
# OccupancyCas_BootDf.names = ["Intercept", "PokeTime","Group:Caspase","PokeTime&Group:Juveniles"]
# OccupancyCas_BootDf.Y = 1:nrow(OccupancyCas_BootDf)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig7","1000CasOccupancyBootstrap.csv"),OccupancyCas_BootDf)
OccupancyCas_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig7","1000CasOccupancyBootstrap.csv"),DataFrame)
OccupancyCas_BootDf.err = [tuple(parse.(Float64,split(err[2:end-1],", "))...) for err in OccupancyCas_BootDf.err]
yprop = ("",font(10, "Bookman Light"),(collect(1:nrow(OccupancyCas_BootDf)),
    [L"Intercept",L"PokeTime",L"Group:Caspase",L"PokeTime\&Group:Caspase"]))
    xprop = ("Coefficient estimate", font(10, "Bookman Light"), (-1.1,0.2))
    @df OccupancyCas_BootDf scatter(:coef ,1:nrow(OccupancyCas_BootDf), xerror = :err,legend = false,
    xaxis = xprop, yaxis = yprop, markercolor = :gray75, size = (600,600))
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
    annotate!([(-1,1,"n.s", 10), (-1,2,"*", 10), (-1,3,"n.s.", 10), (-1,4,"*", 10)])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "J","VirusOccupancyBootstrap.pdf"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","Fig3", "J","VirusOccupancyBootstrap.csv"), select(OccupancyCas_BootDf, Not([:type,:group,:Y])))
########################### SFig 3 ###########################
########################### A ###########################
## Scatter Poke & Trial
@df Age_p scatter(:LogOut, :Streak, markersize = 1, xlabel = "Poke Time(Log10 seconds)",
    ylabel = "Trial", legend = false)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","SFig3", "A","ScatterPoke&Trial.pdf"))
prov = select(Age_p,[:Day,:MouseID,:Age, :Poke, :Streak, :LogOut])
rename(prov, :Streak => :Trial)
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","SFig3", "A","ScatterPoke&Trial.csv"), prov)
########################### B ###########################
## Median Leaving Time
Age_TrialMedianSurvival=plot();
dd1 = combine(groupby(Age_s,[:MouseID,:LongBinnedStreak,:Age]), :LogDuration => median => :LogDuration)
dd2 = combine(groupby(dd1,:Age)) do dd
    group_summary(dd,:LongBinnedStreak,:LogDuration; normality = false)
end
filter!(r-> r.LongBinnedStreak <= 90, dd2)
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","SFig3", "B","MedianRatePerTrial.csv"), dd2)
xprop = ("Trial", (10,90), xyfont)
yprop = ("Median laeving time(seconds)", xyfont,(log(1),log10(30)),(log10.([1,5,10,15,20, 30]),string.([1,5,10,15,20,30])))
@df dd2 scatter(:LongBinnedStreak, :Central, yerror = :ERR, group = :Age, xaxis = xprop, yaxis = yprop, grid = true,size = (600,600), legend = true)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","SFig3", "B","MedianRatePerTrial.pdf"))
########################### C ###########################
##  Leving rate ######################################
Sex_LevingAn, Sex_LevingAn_df = function_analysis(Age_s,:LogDuration, cumulative_algorythm;
    grouping = :Sex, calc = :bootstrapping, color = [:violetred3 :sandybrown])
    xprop = ("Poke Time(seconds)", xyfont,(log10.([0.1,1,10,100,1000]),["0.1","1","10","100","1000"]))
    yprop = ("Probablity of leaving", xyfont)
    plot!(Sex_LevingAn, xaxis = xprop, yaxis = yprop, legend = (0.8, 0.25))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","SFig3", "C","SexLevingRate.pdf"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","SFig3", "C","SexLevingRate.csv"), Sex_LevingAn_df)
########################### D ###########################
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
    xaxis = xprop, yaxis = yprop, markercolor = :gray75)
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","SFig3", "D","SexLevingRateBootstrap.pdf"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","SFig3", "D","SexLevingRateBootstrap.csv"), Sex_btdf)
