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
##### Selection criteria Caspase ######
# AL over trials
Cas_s[!,:BinnedStreak] = bin_axis(Cas_s.Streak; unit_step = 4)
res3 = summary_xy(Cas_s,:BinnedStreak,:AfterLast; group = :Virus)
cas_alXtrial = @df res3 plot(string.(:BinnedStreak),:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, xlabel = "Trial", ylabel = "Pokes after last reward")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","alXtrials.png"))

# trials over time
Cas_s[!,:BinnedStart] = bin_axis(Cas_s.Start./60; unit_step = 2)
res3 = summary_xy(Cas_s,:BinnedStart,:Streak; group = :Virus)
cas_trialXtime = @df res3 plot(string.(:BinnedStart),:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, legend = :topleft, xlabel = "Time (min)",
    ylabel = "Number of trials")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","trialsXtime.png"))

# AL distribution
df1 = group_kde(Cas_s,:AfterLast; group = [:Virus], points = 50)
cas_AlDist = @df df1 plot(:Xaxis,:Mean, ribbon = :Sem, xlims = (0,25),
    linewidth = 1, linecolor = :auto, group = :Virus,
    xlabel = "Pokes after last reward", ylabel = "PDF")
q95 = quantile(collect(skipmissing(Cas_s.AfterLast)),0.95)
vline!([q95],line = :dash, label = "95th percentile")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","pdfXal.png"))

# Interpoke distribution
limit = quantile(collect(skipmissing(Cas_p.PreInterpoke)),0.95)
df1 = group_kde(Cas_p,:PreInterpoke; group = :Virus, points = 1000)
cas_IPDist = @df df1 plot(:Xaxis,:Mean, ribbon = :Sem, linecolor = :auto, xlims = (0,60),
    xlabel = "Inter-poke interval (s)", ylabel = "PDF", group = :Virus)
vline!([limit], line = :dash, label = "95th quantile")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","InterpokePDF.png"))

# Travel distribution
limit = quantile(collect(skipmissing(Cas_s.Travel_to)),0.95)
df1 = group_kde(Cas_s,:Travel_to; points = 1000, group = :Virus)
cas_TrDist = @df df1 plot(:Xaxis,:Mean, ribbon = :Sem, linecolor = :auto, xlims = (0,120),
    xlabel = "Travel time (s)", ylabel = "PDF", group = :Virus)
vline!([limit], line = :dash, label = "95th quantile")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","Travel_toPDF.png"))

# Trial Duration
x = Cas_s.Trial_duration#[Cas_s.Trial_duration .<= 60]
mix_trial_duration = mixture_gamma(x)
interval = 0:0.2:60
plot(interval,pdf(mix_trial_duration,interval),xlims = (0,60), label = "Convex combination",
    yrotation = 60)
plt = twinx()
histogram!(plt,x[x .<=60], nbins = 150, color = :grey, fillalpha = 0.3, xlims = (0,60), linewidth = 0, label = false,
    xlabel = "Trial duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","travelHist_PDF.png"))

plot(interval,pdf(mix_trial_duration,interval),xlims = (0,60), label = "Convex combination")
plot!(interval,pdf(mix_trial_duration.components[1],interval), xlims = (0,60), linecolor = :cyan, label = "First component")
plot!(interval,pdf(mix_trial_duration.components[2],interval), xlims = (0,60),linecolor = :magenta, label = "Second component")
vline!([30], linestyle = :dot, label = "Threshold")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","MixedTravelDists.png"))
## Afterlast df selection
limit_cas = quantile(collect(skipmissing(Cas_s.AfterLast)),0.95)
cas_df = filter(r->
    r.Trial_duration < 30 &&
    r.AfterLast < limit_cas &&
    r.Gen == "Rbp4-cre"
    ,Cas_s)
## AfterLast with Caspase
cas_afterlast = DoubleAnalysis(cas_df,:Virus,:AfterLast, yspan = (0,3))
cas_afterlast.JarqueBera
cas_afterlast.nonparametric_plot
ylabel!("Pokes after last reward")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","ALcaspase.png"))
# cas_afterlast.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + (1|MouseID)),cas_df))
fm2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Virus + (1|MouseID)),cas_df))
Likelyhood_Ratio_test(fm1,fm2)
cas_afterlast.mice_summary
cas_afterlast.parametric_summary
#=
Likelihood ratio test on
AfterLast ~ 1 + Virus + (1|MouseID) versus AfterLast ~ 1 + (1|MouseID): p < 0.01,
effect of Caspase in number of attempts after last reward = -1.26 ± 0.39
=#
## Incorrect with Caspase
cas_correct = DoubleAnalysis(cas_df,:Virus,:IncorrectLeave, yspan = (0,1))
cas_correct.JarqueBera
cas_correct.nonparametric_plot
ylabel!("Fraction of error trials")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","Errorcaspase.png"))
fm1 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + (1|MouseID)),cas_df))
fm2 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + Virus + (1|MouseID)),cas_df))
Likelyhood_Ratio_test(fm1,fm2)
cas_correct.parametric_summary
#=
Likelihood ratio test on
IncorrectLeave ~ 1 + Age + (1|MouseID) versus IncorrectLeave ~ 1 + (1|MouseID): p < 0.01,
effect of Juveniles on probability of correct leave = -0.12 ± 0.03
=#
## Num Rewards with Caspase
cas_rewards = DoubleAnalysis(cas_df,:Virus,:Num_Rewards,summary_opt = :SUM)
cas_rewards.JarqueBera
cas_rewards.nonparametric_plot
ylabel!("Total Rewards")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","Rewardscaspase.png"))
## Interpoke with Caspase
limit = quantile(collect(skipmissing(Cas_p.PreInterpoke)),0.95)
f_cas = filter(r ->
    r.PreInterpoke < limit
    ,Cas_p)
cas_interpoke = DoubleAnalysis(f_cas,:Virus,:PreInterpoke)
cas_interpoke.JarqueBera
cas_interpoke.nonparametric_plot
ylabel!("Inter-poke interval(s)")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","InterpokeCaspase.png"))
fm1 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + (1|MouseID)),f_cas))
fm2 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Virus + (1|MouseID)),f_cas))
fm3 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Virus + Sex + (1|MouseID)),f_cas))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
cas_interpoke.parametric_summary
#=
Likelihood ratio test on
Interpoke ~ 1 + Virus + (1|MouseID) versus Interpoke ~ 1 + (1|MouseID): n.s.
=#
#### Fig3 ###
plot(cas_afterlast.nonparametric_plot,
    cas_correct.nonparametric_plot,
    cas_rewards.nonparametric_plot,
    cas_interpoke.nonparametric_plot,
    thickness_scaling = 1)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","Fig3.png"))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","Fig3.pdf"))
## Travel_to with Caspase
# @df f_cas density(:Travel_to)
# limit = median(Cas_s.Travel_to)
limit = quantile(collect(skipmissing(Cas_s.Travel_to)),0.95)
f_cas = filter(r ->
    r.Travel_to < limit
    ,Cas_s)
cas_travel = DoubleAnalysis(f_cas,:Virus,:Travel_to)
cas_travel.JarqueBera
cas_travel.nonparametric_plot
ylabel!("Travel time (s)")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","TravelCaspase.png"))

cas_travel.parametric_plot
f_cas[!,:Travel_to] = map(Float64,f_cas.Travel_to)
fm1 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + (1|MouseID)),f_cas))
fm2 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Virus + (1|MouseID)),f_cas))
fm3 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Virus + Sex + (1|MouseID)),f_cas))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
cas_travel.parametric_summary
#=
Likelihood ratio test on
Travel_to ~ 1 + Virus + (1|MouseID) versus Interpoke ~ 1 + (1|MouseID): n.s.
=#
###
plot(cas_alXtrial,cas_trialXtime,
    cas_AlDist,cas_IPDist,
    cas_TrDist,cas_travel.nonparametric_plot,
    layout = grid(3,2),
    thickness_scaling = 1,
    size=(600,900))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","SFig3.png"))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig3","SFig3.pdf"))



## Number of Reward
age_rew = DoubleAnalysis(age_df,:Age,:Num_Rewards, summary_opt = :SUM)
age_rew.JarqueBera
age_rew.nonparametric_plot
age_rew.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + (1|MouseID)),age_df))
fm2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + (1|MouseID)),age_df))
fm3 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + Sex + (1|MouseID)),age_df))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
age_rew.parametric_summary
res = age_rew.mice_summary
res[!,:Sex] = [x in females ? "F" : "M" for x in res.MouseID]
res[!,:Combo] = res.Sex .* res.Age
group_summary(res,:Combo, :AfterLast)
open_html_table(age_df)
