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
##### Selection criteria Jeveniles ######
# AL over trials
Age_s[!,:BinnedStreak] = bin_axis(Age_s.Streak; unit_step = 4)
res3 = summary_xy(Age_s,:BinnedStreak,:AfterLast; group = :Age)
age_alXtrial = @df res3 plot(string.(:BinnedStreak),:Mean, group = :Age, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, xlabel = "Trial", ylabel = "Pokes after last reward")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","alXtrials.png"))

# trials over time
Age_s[!,:BinnedStart] = bin_axis(Age_s.Start./60; unit_step = 2)
res3 = summary_xy(Age_s,:BinnedStart,:Streak; group = :Age)
age_trialXtime = @df res3 plot(string.(:BinnedStart),:Mean, group = :Age, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, legend = :topleft, xlabel = "Time (min)",
    ylabel = "Number of trials")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","trialsXtime.png"))

# AL distribution
df1 = group_kde(Age_s,:AfterLast; group = [:Age], points = 50)
age_AlDist = @df df1 plot(:Xaxis,:Mean, ribbon = :Sem, xlims = (0,25),
    linewidth = 1, linecolor = :auto, group = :Age,
    xlabel = "Pokes after last reward", ylabel = "PDF")
q95 = quantile(collect(skipmissing(Age_s.AfterLast)),0.95)
vline!([q95],line = :dash, label = "95th percentile")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","pdfXal.png"))

# Interpoke distribution
limit = quantile(collect(skipmissing(Age_p.PreInterpoke)),0.95)
df1 = group_kde(Age_p,:PreInterpoke; group = :Age, points = 1000)
age_IPDist = @df df1 plot(:Xaxis,:Mean, ribbon = :Sem, linecolor = :auto, xlims = (0,60),
    xlabel = "Inter-poke interval (s)", ylabel = "PDF", group = :Age)
vline!([limit], line = :dash, label = "95th percentile")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","InterpokePDF.png"))

# Travel distribution
limit = quantile(collect(skipmissing(Age_s.Travel_to)),0.95)
df1 = group_kde(Age_s,:Travel_to; points = 1000, group = :Age)
age_TrDist = @df df1 plot(:Xaxis,:Mean, ribbon = :Sem, linecolor = :auto, xlims = (0,120),
    xlabel = "Travel time (s)", ylabel = "PDF", group = :Age)
vline!([limit], line = :dash, label = "95th quantile")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","Travel_toPDF.png"))

###

# Trial Duration
x = Age_s.Trial_duration#[Age_s.Trial_duration .<= 60]
mix_trial_duration = mixture_gamma(x)
interval = 0:0.2:60
plot(interval,pdf(mix_trial_duration,interval),xlims = (0,60), label = "Convex combination",
    yrotation = 60)
plt = twinx()
histogram!(plt,x[x .<=60], nbins = 150, color = :grey, fillalpha = 0.3, xlims = (0,60), linewidth = 0, label = false,
    xlabel = "Trial duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","travelHist_PDF.png"))

plot(interval,pdf(mix_trial_duration,interval),xlims = (0,60), label = "Convex combination")
plot!(interval,pdf(mix_trial_duration.components[1],interval), xlims = (0,60), linecolor = :cyan, label = "First component")
plot!(interval,pdf(mix_trial_duration.components[2],interval), xlims = (0,60),linecolor = :magenta, label = "Second component")
vline!([30], linestyle = :dot, label = "Threshold")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","MixedTravelDists.png"))

## Afterlast df selection
limit_age = quantile(collect(skipmissing(Age_s.AfterLast)),0.95)
age_df = filter(r->
    r.Trial_duration < 30 &&
    r.AfterLast < limit_age
    ,Age_s)
## AfterLast with Juveniles
age_afterlast = DoubleAnalysis(age_df,:Age,:AfterLast)
age_afterlast.JarqueBera
age_afterlast.mice_summary
age_afterlast.nonparametric_plot
ylabel!("Pokes after last reward")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","ALjuvenile.png"))
age_afterlast.nonparametric_summary
# age_afterlast.parametric_plot
# age_afterlast.parametric_summary

fm1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + (1|MouseID)),age_df))
fm2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + (1|MouseID)),age_df))
fm3 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Age + Sex + (1|MouseID)),age_df))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)

fm4 = fit(MixedModel,@formula(AfterLast ~ 1 + (1|MouseID)),age_df,Poisson())
fm5 = fit(MixedModel,@formula(AfterLast ~ 1 + Age + (1|MouseID)),age_df,Poisson())
fm6 = fit(MixedModel,@formula(AfterLast ~ 1 + Age + Sex + (1|MouseID)),age_df,Poisson())
Likelyhood_Ratio_test(fm4,fm5)
Likelyhood_Ratio_test(fm5,fm6)

scatter(age_df.AfterLast,predict(fm5),markersize=2, legend  = false)
p = predict(fm6)
y = age_df.AfterLast
r = @. sign(y - p) * sqrt(devresid(Poisson(), y, p))
scatter(r,age_df.AfterLast;markersize=1, stroke = 1, legend = false)
qqnorm(r; qqline = :R, markersize = 2)
fitted(fm5)
# open_html_table(age_afterlast.mice_summary)
# open_html_table(filter(r-> r.MouseID == "RJ67", age_df))
#=
Likelihood ratio test on
AfterLast ~ 1 + Age + (1|MouseID) versus AfterLast ~ 1 + (1|MouseID): p < 0.01,
effect of Juveniles in number of attempts after last reward = -0.70 ± 0.19
=#
## Incorrect with Juveniles
age_correct = DoubleAnalysis(age_df,:Age,:IncorrectLeave, yspan = (0,1))
age_correct.JarqueBera
age_correct.nonparametric_plot
# age_correct.parametric_plot
ylabel!("Fraction of error trials")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","Errorjuvenile.png"))

fm1 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + (1|MouseID)),age_df))
fm2 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + Age + (1|MouseID)),age_df))
fm3 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + Age + Sex + (1|MouseID)),age_df))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
age_correct.parametric_summary
res = age_correct.mice_summary
res[!,:Sex] = [x in females ? "F" : "M" for x in res.MouseID]
res[!,:Combo] = res.Sex .* res.Age
group_summary(res,:Combo, :IncorrectLeave)
JarqueBeraTest(residuals(fm2))
# open_html_table(age_correct.mice_summary)
#=
Likelihood ratio test on
IncorrectLeave ~ 1 + Age + (1|MouseID) versus IncorrectLeave ~ 1 + (1|MouseID): p < 0.01,
effect of Juveniles on probability of correct leave = -0.07 ± 0.02
=#
## Num Rewards with Juveniles
age_rewards = DoubleAnalysis(age_df,:Age,:Num_Rewards,summary_opt = :SUM)
age_rewards.JarqueBera
age_rewards.nonparametric_plot
ylabel!("Total Rewards")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","Rewardsjuvenile.png"))
## Interpoke with Juveniles
limit = quantile(collect(skipmissing(Age_p.PreInterpoke)),0.95)
f_age = filter(r ->
    r.PreInterpoke < limit
    ,Age_p)
age_interpoke = DoubleAnalysis(f_age,:Age,:PreInterpoke)
age_interpoke.JarqueBera
age_interpoke.nonparametric_plot
ylabel!("Inter-poke interval(s)")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","InterpokeEffect.png"))

age_interpoke.parametric_plot
fm1 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + (1|MouseID)),f_age))
fm2 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Age + (1|MouseID)),f_age))
fm3 = fit!(LinearMixedModel(@formula(PreInterpoke ~ 1 + Age + Sex + (1|MouseID)),f_age))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
age_interpoke.parametric_summary

#### Fig1 ###
plot(age_afterlast.nonparametric_plot,
    age_correct.nonparametric_plot,
    age_rewards.nonparametric_plot,
    age_interpoke.nonparametric_plot,
    thickness_scaling = 1)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","Fig1.png"))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","Fig1.pdf"))

## Travel_to with Juveniles
limit = quantile(collect(skipmissing(Age_s.Travel_to)),0.95)
f_age = filter(r ->
    r.Travel_to < limit
    ,Age_s)
age_travel = DoubleAnalysis(f_age,:Age,:Travel_to)
age_travel.JarqueBera
age_travel.nonparametric_plot
ylabel!("Travel time (s)")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","TravelEffect.png"))

age_travel.parametric_plot
f_age[!,:Travel_to] = map(Float64,f_age.Travel_to)
fm1 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + (1|MouseID)),f_age))
fm2 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Age + (1|MouseID)),f_age))
fm3 = fit!(LinearMixedModel(@formula(Travel_to ~ 1 + Age + Sex + (1|MouseID)),f_age))
Likelyhood_Ratio_test(fm1,fm2)
Likelyhood_Ratio_test(fm2,fm3)
age_travel.parametric_summary
###
plot(age_alXtrial,age_trialXtime,
    age_AlDist,age_IPDist,
    age_TrDist,age_travel.nonparametric_plot,
    layout = grid(3,2),
    thickness_scaling = 1,
    size=(600,900))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","SFig1.png"))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","SFig1","SFig1.pdf"))
