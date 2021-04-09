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
FAge_s[!,:PokingRatio] = FAge_s.Poking_time./FAge_s.Trial_duration
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
FCas_s[!,:PokingRatio] = FCas_s.Poking_time./FCas_s.Trial_duration
RatioCas0 = fit(MixedModel,@formula(PokingRatio ~ 1 + Streak + (1|MouseID)),FCas_s)
RatioCas1 = fit(MixedModel,@formula(PokingRatio ~ 1 + Streak + Virus + (1|MouseID)),FCas_s)
Likelyhood_Ratio_test(RatioCas0,RatioCas1)



##
differenza = Matrix(reshapecas[:, Not(:Bin_Streak)]) - Matrix(reshapetom[:, Not(:Bin_Streak)])
heatmap(xlab, ylab, differenza,
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","5-HeatmapDiff.pdf"))
###############################################################################
###############################################################################
###############################################################################
transform!(Cas_p, [:Streak, :Out] .=> zscore)
verbagg1 = @formula(Leave ~ 1 + Streak_zscore * Virus + Out_zscore * Virus +  (1|MouseID));
mdl1 = fit(MixedModel,verbagg1, Cas_p, Bernoulli())
Cas_p.Model_Leave = predict(mdl1)
rng = MersenneTwister(1234321)
samp1 = parametricbootstrap(rng,100,mdl1)
df = DataFrame(samp1.allpars)
bootdf = combine(groupby(df,[:type, :group, :names]), :value => shortestcovint => :interval)

bootdf.coef = push!(coef(mdl1), mean(ranef(mdl1)[1]))
bootdf.variable = ["Intercept", "Streak", "Virus:Caspase", "Out",
    "Streak & Virus: Caspase",
    "Out & Virus: Caspase", "MouseID"]
bootdf.err = [(x[2],x[1]) for (x,c) in zip(bootdf.interval,boot]
transform!(bootdf, [:coef, :interval] => ByRow((c,e) -> (c -e[1], e[2]-c)) => :err)

open_html_table(bootdf)
@df bootdf[1:end-1,:] scatter(:coef ,1:nrow(bootdf)-1,
    # xerror = :interval,
    xerror = :err,
    yticks = (1:nrow(bootdf), :variable))
vline!([0], linecolor = :red, legend = false)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","10-Model.pdf"))
###############
df0 = filter(r -> r.Streak <= 50, Cas_p)
transform!(df0, :LogOut => (x -> bin_axis(x; length = 6)) => :Bin_LogOut)
transform!(df0, :Streak => (x -> bin_axis(x; unit_step = 10)) => :Bin_Streak)
df1 = combine(groupby(df0,[:MouseID,:Bin_LogOut,:Virus,:Bin_Streak]),
    [:Leave, :Model_Leave] .=> mean .=> [:Leave, :Model_Leave])
df2 = combine(groupby(df1,[:Bin_LogOut,:Virus,:Bin_Streak]),
    [:Leave, :Model_Leave] .=> mean .=> [:Leave, :Model_Leave])
##

dfcas = filter(r -> r.Virus == "Caspase", df2)
filter!(r -> -1 <= r.Bin_LogOut <= 2.1,dfcas)
reshapecas = unstack(dfcas,:Bin_Streak, :Bin_LogOut, :Model_Leave)
dropmissing!(reshapecas)
ylab = string.(sort(union(reshapecas.Bin_Streak)))
xlab = string.(sort(union(dfcas.Bin_LogOut)))
heatmap(xlab, ylab, Matrix(reshapecas[:, Not(:Bin_Streak)]), clim = (0.15,0.38),
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","7-HeatmapCas.pdf"))
##
dfcas = filter(r -> r.Virus == "tdTomato", df2)
filter!(r -> -1 < r.Bin_LogOut <= 2.1,dfcas)
reshapetom = unstack(dfcas,:Bin_Streak, :Bin_LogOut, :Model_Leave)
dropmissing!(reshapetom)
ylab = string.(sort(union(reshapetom.Bin_Streak)))
xlab = string.(sort(union(dfcas.Bin_LogOut)))
heatmap(xlab, ylab, Matrix(reshapetom[:, Not(:Bin_Streak)]),clim = (0.15,0.38),
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","8-HeatmapTomato.pdf"))
##
differenza = Matrix(reshapecas[:, Not(:Bin_Streak)]) - Matrix(reshapetom[:, Not(:Bin_Streak)])
heatmap(xlab, ylab, differenza,
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","9-HeatmapDiff.pdf"))
