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
