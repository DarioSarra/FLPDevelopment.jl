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
##
Age_Rewards = Difference(Age_s, :Age, :Num_Rewards, ylabel = "number of rewards per trial", ylims = (0,1.3))
Age_Rewards.plt
Age_Rewards.test
Age_Rewards.groupdf
Age_Rewards.groupdf[!,:Measure] .= "Number of rewards"
rename!(Age_Rewards.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Post-Rewievs", "Rewards.pdf"
##
NumPokes_f = @formula(Leave ~ 1 + PokeInStreak +  (1+PokeInStreak|MouseID));
NumPokes_m = fit(MixedModel,NumPokes_f, Age_p, Bernoulli())
PokeTime_f = @formula(Leave ~ 1 + LogOut_zscore +  (1+LogOut_zscore|MouseID));
PokeTime_m = fit(MixedModel,PokeTime_f, Age_p, Bernoulli())
both_f = @formula(Leave ~ 1 + LogOut_zscore +  PokeInStreak +(1+LogOut_zscore+PokeInStreak|MouseID))
both_m = fit(MixedModel,both_f, Age_p, Bernoulli())
AIC_test(PokeTime_m, NumPokes_m)
aic(PokeTime_m)
aic(NumPokes_m)
MixedModels.likelihoodratiotest(PokeTime_m, both_m)
MixedModels.likelihoodratiotest(NumPokes_m, both_m)
))
##
Age_AfterLast = Difference(Age_s, :Age, :AfterLast,
    ylabel = "Consecutive failures before leaving", ylims = (0,5))
Age_AfterLast.plt
Age_AfterLast.test
Age_AfterLast.groupdf
Age_AfterLast.groupdf[!,:Measure] .= "ConsecutiveFailures"
rename!(Age_AfterLast.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Post-Rewievs", "ConsecutiveFailures.pdf"))
##
Cas_Rewards = Difference(Cas_s, :Virus, :Num_Rewards, ylabel = "number of rewards per trial", ylims = (0,1.3))
Cas_Rewards.plt
Cas_Rewards.test
Cas_Rewards.groupdf
Cas_Rewards.groupdf[!,:Measure] .= "Number of rewards"
rename!(Cas_Rewards.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Post-Rewievs", "Cas_Rewards.pdf"))
##
Cas_AfterLast = Difference(Cas_s, :Virus, :AfterLast,
    ylabel = "Consecutive failures before leaving", ylims = (0,5))
Cas_AfterLast.plt
Cas_AfterLast.test
Cas_AfterLast.groupdf
Cas_AfterLast.groupdf[!,:Measure] .= "ConsecutiveFailures"
rename!(Cas_AfterLast.groupdf,[:Central => :Median, :ERR => :CI])
HypothesisTests.OneSampleTTest(mousecoeff.Num_Rewards .+ stimbound_m.β[2])
