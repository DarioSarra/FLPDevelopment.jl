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
## Fig 1
age_al_lim = quantile(Age_s.AfterLast,0.95)
age_df = filter(r->
    # r.Streak <=80
    # r.Performance >25 &&
    r.Trial_duration < 30 &&
    r.AfterLast <= age_al_lim &&
    r.Stop/60 < 50
    ,Age_s)
age_afterlast = DoubleAnalysis(age_df,:Age,:AfterLast)
ylabel!(age_afterlast.nonparametric_plot,"Pokes after last reward")
age_correct = DoubleAnalysis(age_df,:Age,:IncorrectLeave, yspan = (0,1))
ylabel!(age_correct.nonparametric_plot,"Fraction of errors")
# age_rewards = DoubleAnalysis(age_df,:Age,:Num_Rewards,summary_opt = :SUM)
# ylabel!(age_rewards.nonparametric_plot,"Number of rewards")
age_rewrate = DoubleAnalysis(age_df,:Age,:RewRate)
ylabel!(age_rewrate.nonparametric_plot,"Reward rate")
age_limit_interpoke = quantile(collect(skipmissing(Age_p.PreInterpoke)),0.95)
age_f_interpoke = filter(r -> 0 < r.PreInterpoke < 1, Age_p)
nrow(age_f_interpoke)/ nrow(Age_p)
age_interpoke = DoubleAnalysis(age_f_interpoke,:Age,:PreInterpoke)
ylabel!(age_interpoke.nonparametric_plot,"Inter poke interval (s)")
## Fig 3
cas_al_lim = quantile(Cas_s.AfterLast,0.95)
cas_df = filter(r->
    # r.Streak <=80 &&
    # r.Performance >25 &&
    r.Trial_duration < 30 &&
    r.AfterLast < cas_al_lim &&
    r.Stop/60 <= 50  &&
    r.Gen == "Rbp4-cre"
    ,Cas_s)
cas_afterlast = DoubleAnalysis(cas_df,:Virus,:AfterLast)
ylabel!(cas_afterlast.nonparametric_plot,"Pokes after last reward")
cas_correct = DoubleAnalysis(cas_df,:Virus,:IncorrectLeave, yspan = (0,1))
ylabel!(cas_correct.nonparametric_plot,"Fraction of errors")
# cas_rewards = DoubleAnalysis(cas_df,:Virus,:Num_Rewards,summary_opt = :SUM)
# ylabel!(cas_rewards.nonparametric_plot,"Number of rewards")
cas_rewrate = DoubleAnalysis(cas_df,:Virus,:RewRate)
ylabel!(cas_rewrate.nonparametric_plot,"Reward rate")
cas_f_interpoke = filter(r -> 0 < r.PreInterpoke < 1, Cas_p)
cas_interpoke = DoubleAnalysis(cas_f_interpoke,:Virus,:PreInterpoke)
ylabel!(cas_interpoke.nonparametric_plot,"Inter poke interval (s)")
## Fig S1
# AL distribution
df1 = group_cdf(Age_s,:AfterLast)
age_AlCum = @df df1 plot(:Xaxis,:Mean, ribbon = :Sem, xlims = (0,25),
    linewidth = 1, linecolor = :auto, legend = false, xticks = 0:2:26,
    xlabel = "Pokes after last reward", ylabel = "Cumulative distribution")
vline!([age_al_lim],line = :dash)
hline!([0.95], line=:dash, label = "95th percentile")
# trials over time
Age_s[!,:BinnedStart] = bin_axis(Age_s.Start./60; unit_step = 2)
res3 = summary_xy(Age_s,:BinnedStart,:Streak; group = :Age)
age_trialXtime = @df res3 plot(string.(:BinnedStart),:Mean, group = :Age, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, legend = :topleft, xlabel = "Time (min)",
    ylabel = "Number of trials")
vline!([25], line = :dash, label = "Threshold")
# AL over trials
age_df[!,:BinnedStreak] = bin_axis(age_df.Streak; unit_step = 4)
res3 = summary_xy(age_df,:BinnedStreak,:AfterLast; group = :Age)
age_alXtrial = @df res3 plot(string.(:BinnedStreak),:Mean, group = :Age, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, xlabel = "Trial", ylabel = "Pokes after last reward")
# cumulative reward
res3 = summary_xy(age_df,:Streak,:Cum_Rewards; group = :Age)
filter!(r -> r.Streak <=65 ,res3)
age_rewXtrials = @df res3 plot(:Streak,:Mean, group = :Age, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, legend = :topleft, xlabel = "Trials",
    ylabel = "Reward")
# travel time
age_limit_travel = quantile(age_df.Travel_to[age_df.PreInterpoke.> 0],0.95)
age_f_travel = filter(r -> 0 < r.Travel_to < age_limit_travel, age_df)
age_travel = DoubleAnalysis(age_f_travel,:Age,:Travel_to)
ylabel!(age_travel.nonparametric_plot,"Travel time (s)")
## Fig S3
df1 = group_cdf(Cas_s,:AfterLast)
cas_AlCum = @df df1 plot(:Xaxis,:Mean, ribbon = :Sem, xlims = (0,25),
    linewidth = 1, linecolor = :auto, legend = false, xticks = 0:2:26,
    xlabel = "Pokes after last reward", ylabel = "Cumulative distribution")
vline!([cas_al_lim],line = :dash)
hline!([0.95], line=:dash, label = "95th percentile")
# trials over time
Cas_s[!,:BinnedStart] = bin_axis(Cas_s.Start./60; unit_step = 2)
res3 = summary_xy(Cas_s,:BinnedStart,:Streak; group = :Virus)
cas_trialXtime = @df res3 plot(string.(:BinnedStart),:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, legend = :topleft, xlabel = "Time (min)",
    ylabel = "Number of trials")
vline!([25], line = :dash, label = "Threshold")
# AL over trials
cas_df[!,:BinnedStreak] = bin_axis(cas_df.Streak; unit_step = 4)
res3 = summary_xy(cas_df,:BinnedStreak,:AfterLast; group = :Virus)
cas_alXtrial = @df res3 plot(string.(:BinnedStreak),:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, xlabel = "Trial", ylabel = "Pokes after last reward")
# travel time
cas_travel = DoubleAnalysis(cas_df,:Virus,:Travel_to)
ylabel!(cas_travel.nonparametric_plot,"Travel time (s)")
cas_travel = DoubleAnalysis(cas_df,:Virus,:Travel_to, summary_opt = :MEDIAN)
ylabel!(cas_travel.nonparametric_plot,"Travel time (s)")
# cumulative reward
res3 = summary_xy(cas_df,:Streak,:Cum_Rewards; group = :Virus)
filter!(r -> r.Streak <=65 ,res3)
cas_rewXtrials = @df res3 plot(:Streak,:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, legend = :topleft, xlabel = "Trials",
    ylabel = "Reward")
##
Fig1 = plot(age_afterlast.nonparametric_plot,
    age_correct.nonparametric_plot,
    age_rewrate.nonparametric_plot,
    # age_rewards.nonparametric_plot,
    age_interpoke.nonparametric_plot,
    # age_travel.nonparametric_plot,
    thickness_scaling = 1)
##
savefig(Fig1,joinpath(replace(path,basename(path)=>""),"Development_Figures","ShortPlotting","Fig1.pdf"))
##
Fig3 = plot(cas_afterlast.nonparametric_plot,
    cas_correct.nonparametric_plot,
    # cas_rewards.nonparametric_plot,
    age_rewrate.nonparametric_plot,
    cas_interpoke.nonparametric_plot,
    # cas_travel.nonparametric_plot,
    thickness_scaling = 1)
##
savefig(Fig3,joinpath(replace(path,basename(path)=>""),"Development_Figures","ShortPlotting","Fig3.pdf"))
##
FigS1 = plot(age_AlCum,
    age_trialXtime,
    # age_alXtrial,
    age_rewXtrials,
    # age_travel.nonparametric_plot,
    age_rewards.nonparametric_plot,
    layout = grid(2,2),
    thickness_scaling = 1,
    size=(600,600))
##
savefig(FigS1, joinpath(replace(path,basename(path)=>""),"Development_Figures","ShortPlotting","FigS1.pdf"))
##
FigS3 = plot(cas_AlCum,
    cas_trialXtime,
    # cas_alXtrial,
    # cas_RewCum,
    cas_rewXtrials,
    # cas_travel.nonparametric_plot,
    cas_rewards.nonparametric_plot,
    layout = grid(2,2),
    thickness_scaling = 1,
    size=(600,600))
##
savefig(FigS3, joinpath(replace(path,basename(path)=>""),"Development_Figures","ShortPlotting","FigS3.pdf"))
