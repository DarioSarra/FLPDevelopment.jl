using Revise, FLPDevelopment, BrowseTables, Random
gr(size=(600,600), tick_orientation = :out, grid = false,
    linecolor = :black,
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
#= to do
    1 - remove pokes after an interpoke longer than 7s or 1s
    2 - collapse pokes with an interpoke shorter than 50 ms (see report zach new protocols: typical span 0.05 to 0.5 seconds)
        find the cases i
        remove :Stim, :StimFreq, :Wall, :Block, :Delta, :Turn, :StreakInBlock,

        check i and i+1 are in the same :Session, :Bout, :Streak, :Side, :Day, :MouseID, :Gen,
            :ExpSession, :ProtocolSession, :Age, :Sex, :Performance,

        from case i take :Poke, :SideHigh, Correct, :PokeHierarchy, :PokeIn, :PreInterpoke, :PokeNum,  :In

        from case i+1 take :Leave, :PokeOut, :PostInterpoke, :Out

        from case i and i+1 take any(, Reward)

        immediately recalculate PokeDur, :PokeInStreak

        recalculate in the end :Poke
=#
Remove_vals = [:Stim, :StimFreq, :Wall, :Block, :Delta, :Turn, :StreakInBlock, :ReverseStreak]
Same_vals = [:Session, :Streak, :Side, :Day, :MouseID, :Gen,
        :ExpSession, :ProtocolSession, :Age, :Sex, :Performance]
PreviousPoke_vals = [:Poke, :SideHigh, :PokeHierarchy, :PokeIn, :PreInterpoke, :PokeNum, :PokeInStreak, :In]
NextPoke_vals = [:Bout, :Leave, :PokeOut, :PostInterpoke, :Out]
# remove useless columns
NAge_p = Age_p[:,Not(Remove_vals)]
FAge_p = combine(groupby(NAge_p, :Session)) do dd
    rm_interpokes(copy(dd); shortlimit = 0.1, longlimit = 50000)
end
FAge_s = reprocess_streaks(FAge_p)
FAge_s.Age = categorical(FAge_s.Age)
levels!(FAge_s.Age, ["Adults", "Juveniles"])
AlAge0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),FAge_s)
AlAge1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Age + (1|MouseID)),FAge_s)
Likelyhood_Ratio_test(AlAge0,AlAge1)
PokingAge0 = fit(MixedModel,@formula(Poking_time ~ 1 + Streak + (1|MouseID)),FAge_s)
PokingAge1 = fit(MixedModel,@formula(Poking_time ~ 1 + Streak + Age + (1|MouseID)),FAge_s)
Likelyhood_Ratio_test(PokingAge0,PokingAge1)
FAge_s[!,:PokingRatio] = FAge_s.Poking_time./FAge_s.Trial_Duration
RatioAge0 = fit(MixedModel,@formula(PokingRatio ~ 1 + Streak + (1|MouseID)),FAge_s)
RatioAge1 = fit(MixedModel,@formula(PokingRatio ~ 1 + Streak + Age + (1|MouseID)),FAge_s)
Likelyhood_Ratio_test(RatioAge0,RatioAge1)
##
NCas_p = Cas_p[:,Not(Remove_vals)]
FCas_p = combine(groupby(NCas_p, :Session)) do dd
    rm_interpokes(copy(dd); shortlimit = 0.1, longlimit = 50000)
end
FCas_s = reprocess_streaks(FCas_p)
FCas_s.Virus = categorical(FCas_s.Virus)
levels!(FCas_s.Virus, ["tdTomato", "Caspase"])
AlCas0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),FCas_s)
AlCas1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),FCas_s)
Likelyhood_Ratio_test(AlCas0,AlCas1)
PokingCas0 = fit(MixedModel,@formula(Poking_time ~ 1 + Streak + (1|MouseID)),FCas_s)
PokingCas1 = fit(MixedModel,@formula(Poking_time ~ 1 + Streak + Virus + (1|MouseID)),FCas_s)
Likelyhood_Ratio_test(PokingCas0,PokingCas1)
FCas_s[!,:PokingRatio] = FCas_s.Poking_time./FCas_s.Trial_Duration
RatioCas0 = fit(MixedModel,@formula(PokingRatio ~ 1 + Streak + (1|MouseID)),FCas_s)
RatioCas1 = fit(MixedModel,@formula(PokingRatio ~ 1 + Streak + Virus + (1|MouseID)),FCas_s)
Likelyhood_Ratio_test(RatioCas0,RatioCas1)
