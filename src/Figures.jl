using Revise, FLPDevelopment, BrowseTables
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
    # r.MouseID != "RJ67" && # biting, see B3_RJ67_2020-09-28 minute 7:33
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
agedf = filter(r -> r.Limit, Age_s)
casdf = filter(r -> r.Limit, Cas_s)
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
agedf[!,:BinnedStreak] = bin_axis(agedf.Streak; unit_step = 5)
casdf[!,:BinnedStreak] = bin_axis(casdf.Streak; unit_step = 5)
nrow(agedf)/nrow(Age_s)
nrow(agedf)-nrow(Age_s)
nrow(casdf)/nrow(Cas_s)
nrow(casdf)-nrow(Cas_s)
open_html_table(FLPDevelopment.summarydf(Age_s,Age_p))
open_html_table(FLPDevelopment.summarydf(Cas_s,Cas_p))
########################### Example Session Plots ######################################
tt = filter(r -> r.MouseID == "RJ23" && 08<= r.Streak <= 52, Age_p)
plt = plot(legend = false, xlims = (-1,16), xlabel = "Trials",
    yaxis = false, yticks = false, ylabel = "Time (seconds)")
for r in eachrow(tt)
    FLPDevelopment.session_plot!(plt, r)
end
plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","Session.pdf"))
########################### Afterlast Plots ######################################
#calculate mean of var for each mouse over a trial's bin
overtrialplot = FLPDevelopment.overtrial_plot(agedf, :Age, :AfterLast)
ylabel!("Mean consecutive failures")
xlabel!("Trials")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","AgeALTrials.pdf"))

#calculate mean frequency of afterlas before leaving per mouse and group
frequencyplot = FLPDevelopment.frequency_plot(agedf, :Age, :AfterLast)
ylabel!("Frequency")
xlabel!("Consecutive failures")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","AgeALFreq.pdf"))

#calculate mean and sem afterlast per mouse and group
FLPDevelopment.mean_sem_scatter(agedf, :Age, :AfterLast)[1]
agedf[:,[:Age,:Sex]]
ylabel!("Mean consecutive failures")
xlabel!("Group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","AgeALScat.pdf"))

#calculate mean cumulative of afterlas before leaving per mouse and group
tt = FLPDevelopment.individual_cdf(agedf,:AfterLast)
cdfplot = FLPDevelopment.cdf_plot(agedf, :Age, :AfterLast; estimated = true)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","AgeALCum.pdf"))

agedf.Age = categorical(agedf.Age)
ALAge0 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak  + (1|MouseID)),agedf))
ALAge1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Age + (1|MouseID)),agedf))
ALAge2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Age + Sex + (1|MouseID)),agedf))
Likelyhood_Ratio_test(ALAge0,ALAge1)
Likelyhood_Ratio_test(ALAge1,ALAge2)
coef(ALAge1)
PALAge0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),agedf,Poisson())
PALAge1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Age + (1|MouseID)),agedf,Poisson())
Likelyhood_Ratio_test(PALAge0,PALAge1)
mm = ALAge1
show(
    DataFrame(Formula = mm.formula, modeldof = dof(mm), deviance = deviance(mm))
)
ccdf(Distributions.Chisq(3),8)
########################### Correct Plots ######################################
FLPDevelopment.incorrect_fraction_scatter(agedf,:Age,:IncorrectLeave)[1]
ylabel!("Mean probability of errors")
xlabel!("Group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","AgeERRScat.pdf"))

check = filter(r -> r.CorrectStart, agedf)
ILAge0 = fit(MixedModel,@formula(CorrectLeave ~ 1 + Streak + (1|MouseID)),check,Bernoulli())
ILAge1 = fit(MixedModel,@formula(CorrectLeave ~ 1 + Streak + Age + (1|MouseID)),check,Bernoulli())
ILAge2 = fit(MixedModel,@formula(CorrectLeave ~ 1 + Streak + Age + Sex + (1|MouseID)),check,Bernoulli())
Likelyhood_Ratio_test(ILAge0,ILAge1)
Likelyhood_Ratio_test(ILAge1,ILAge2)
########################### Afterlast Plots ######################################
#calculate mean of var for each mouse over a trial's bin
overtrialplot = FLPDevelopment.overtrial_plot(casdf, :Virus, :AfterLast)
ylabel!("Mean consecutive failures")
xlabel!("Trials")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","VirusALTrials.pdf"))

#calculate mean frequency of afterlas before leaving per mouse and group
frequencyplot = FLPDevelopment.frequency_plot(casdf, :Virus, :AfterLast)
ylabel!("Frequency")
xlabel!("Consecutive failures")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","VirusALFreq.pdf"))

#calculate mean and sem afterlast per mouse and group
unique(casdf[:,:MouseID])
FLPDevelopment.mean_sem_scatter(casdf, :Virus, :AfterLast)[1]
ylabel!("Mean consecutive failures")
xlabel!("Group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","VirusALScat.pdf"))

#calculate mean cumulative of afterlas before leaving per mouse and group
cdfplot = FLPDevelopment.cdf_plot(casdf, :Virus, :AfterLast; estimated = true)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","VirusALCum.pdf"))

check = filter(r -> r.CorrectStart,Cas_s)
ALCas0 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak  + (1|MouseID)),casdf))
ALCas1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),casdf))
Likelyhood_Ratio_test(ALCas0,ALCas1)
PALCas0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),casdf,Poisson())
PALCas1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),casdf,Poisson())
Likelyhood_Ratio_test(PALCas0,PALCas1)
wd = countmap(Cas_s.AfterLast)
get(wd,5,0)
Cas_s.Weights = [get(wd,x,0) for x in Cas_s.AfterLast]
WALCas0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),Cas_s)
WALCas1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),Cas_s)
fm = @formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID))
fm1 = fit(MixedModel, fm, Cas_s)
StatsBase.leverage(fm1)

Cas0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),Cas_s)
mod = Cas0
tune = 	4.685
resid = residuals(mod)
m = median(resid)
MAD = @. median(abs(resid - m))
s = MAD/0.6745
h = MixedModels.leverage(mod)
r = @. resid/(tune*s*sqrt(1-h))
w = [(abs(rr)<1) ? (1 - rr.^2).^2 : 0.0 for rr in r]

WALCas1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),Cas_s, wts = w)
WALCas1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),Cas_s)



WALCas1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),Cas_s)
Likelyhood_Ratio_test(WALCas0,WALCas1)

########################### Correct Plots ######################################
FLPDevelopment.incorrect_fraction_scatter(casdf,:Virus,:IncorrectLeave)[2]
ylabel!("Mean probability of errors")
xlabel!("Group")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig3","VirusERRScat.pdf"))

check = filter(r -> r.CorrectStart, casdf)
ILCas0 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + Streak  + (1|MouseID)),check))
ILCas1 = fit!(LinearMixedModel(@formula(IncorrectLeave ~ 1 + Streak + Virus + (1|MouseID)),check))
Likelyhood_Ratio_test(ILCas0,ILCas1)
##
cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast; nonparametric = true, yspan = (0,6))
# cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast)
cas_afterlast.plot
cas_afterlast.normality
age_afterlast = dvAnalysis(age_df,:Age,:AfterLast; nonparametric = true)#; yspan = (1,3)
# age_afterlast = dvAnalysis(age_df,:Age,:AfterLast)
age_afterlast.plot
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/AfterLastJuv.pdf")
###################### Probability Plots #########################################
cas_correct = dvAnalysis(cas_df,:Virus,:CorrectLeave,yspan = (0,1))
cas_correct.plot
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/CorrectCas.pdf")
age_correct = dvAnalysis(age_df,:Age,:CorrectLeave, yspan = (0,1))
age_correct.plot
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/CorrectJuv.pdf")

############################ Pokes  example #################################
cases = findall((Cas_p.Reward .== false) .& (Cas_p.Correct .== true) .& (Cas_p.PokeInStreak .== 2))
cases = findall((Age_p.Reward .== false) .& (Age_p.Correct .== true) .& (Age_p.PokeInStreak .== 2))
##
idx = cases[8]
p = plot(;legend = false)
# from row 13 to 19
for i in idx - 3:idx + 3
    poke_plot!(p,Age_p[i,:])
end
p

##################### Probability mass function ####################################
# Caspase
# cas_df = filter(r->
#     r.Gen == "Rbp4-cre"&&
#     r.ProtocolSession == 1
#     ,Cas_s)
xvar = :AfterLast
yvar = :AfterLast_frequency
gdc = groupby(cas_df,[xvar,:MouseID,:Virus])
df1 = combine(yvar => mean => yvar,gdc)
sort!(df1,xvar)
df2 = DataFrame(AfterLast = Int[], Count = Int64[])
for (k,i) in countmap(df1[:,xvar])
    push!(df2, [k,i])
end
sort!(df2,xvar)
limit = maximum(df2[df2.Count .>= floor(maximum(df2.Count) * 0.95),xvar])
df3 = combine(groupby(df1,[xvar,:Virus])) do dd
    m = mean(dd[:,yvar])
    SEM = sem(dd[:,yvar])
    (Central = m, ERR = SEM)
end
df4 = filter(r -> r.AfterLast <= 15 #=limit=#, df3)
df4
@df df4 scatter(:AfterLast,:Central,
    group = :Virus,
    yerror = :ERR,
    yticks = 0:0.05:0.4, xticks = 0:1:60, grid = true,
    linecolor = :auto,
    markersize = 3, legend = false, color_palette = [:black,:red,])
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/SplitPMFCas.pdf")
## Probability mass function Age
age_df = filter(r->
    r.ProtocolSession == 1
    ,Age_s)
xvar = :AfterLast
yvar = :AfterLast_frequency
gdc = groupby(age_df,[xvar,:MouseID,:Age])
df1 = combine(yvar => mean => yvar,gdc)
sort!(df1,xvar)
df2 = DataFrame(AfterLast = Int[], Count = Int64[])
for (k,i) in countmap(df1[:,xvar])
    push!(df2, [k,i])
end
sort!(df2,xvar)
# limit = maximum(df2[df2.Count .>= floor(maximum(df2.Count) * 0.95),xvar])
df3 = combine(groupby(df1,[xvar,:Age])) do dd
    m = mean(dd[:,yvar])
    SEM = sem(dd[:,yvar])
    (Central = m, ERR = SEM)
end
df4 = filter(r -> r.AfterLast <= limit, df3)
df4
@df df4 scatter(:AfterLast,:Central,
    group = :Age,
    yerror = :ERR,
    yticks = 0:0.05:0.4, xticks = 0:1:60, grid = true,
    linecolor = :auto,
    markersize = 3, legend = false, color_palette = [:black,:red])
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/SplitPMFJuv.pdf")
###################### Interpoke interval time ####################################
# Caspase
limit = quantile(collect(skipmissing(Cas_p.PreInterpoke)),0.95)
cas_df = filter(r ->
    r.Gen == "HET" &&
    !ismissing(r.PreInterpoke) &&
    r.PreInterpoke < 1
    ,Cas_p)
cas_interpoke = dvAnalysis(cas_df,:Virus,:PreInterpoke; yspan = (0.1,0.5))
cas_interpoke.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/InterpokeCas.pdf")
## Interpoke interval time Age
limit = quantile(collect(skipmissing(Age_p.PreInterpoke)),0.95)
age_df = filter(r ->
    !ismissing(r.PreInterpoke) &&
    r.PreInterpoke < 1
    ,Age_p)
age_interpoke = dvAnalysis(age_df,:Age,:PreInterpoke; yspan = (0.1,0.5))
age_interpoke.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/InterpokeJuv.pdf")
######################### Simple Travel time #######################################
# Caspase
limit = quantile(collect(skipmissing(Cas_s.Travel_to)),0.95)
cas_df = filter(r ->
    !ismissing(r.Travel_to) &&
    r.Travel_to < limit
    ,Cas_s)
cas_interpoke = dvAnalysis(cas_df,:Virus,:Travel_to; yspan = (8,33))
cas_interpoke.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/SimpleTravelCas.pdf")
## Age
limit = quantile(collect(skipmissing(Age_s.Travel_to)),0.95)
age_df = filter(r ->
    !ismissing(r.Travel_to) &&
    r.Travel_to < limit
    ,Age_s)
age_interpoke = dvAnalysis(age_df,:Age,:Travel_to; yspan = (8,33))
age_interpoke.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/SimpleTravelJuv.pdf")
######################### Complex Travel time #######################################
#Caspase
pre_cas = filter(r -> 0 < r.Travel_to <= 60, Cas_s)
pre_cas[!,:Travel_to] = round.(pre_cas.Travel_to, digits = 2)
x = pre_cas.Travel_to
mix_travel = mixture_gamma(x)
p = @df pre_cas histogram(:Travel_to, nbins = 100,label = "data histogram")
m = p.series_list[1].plotattributes[:y]
scaling = maximum(filter(!isnan,m))/maximum(pdf(mix_travel,0:60))
plot!(0:60,pdf(mix_travel,0:60)*scaling, label = "mixture distribution fit", linecolour = :red)
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/Travel_Mix_Gamma_fit_Cas.pdf")
##
@df pre_cas histogram(:Travel_to, nbins = 100, color = :grey, fillalpha = 0.3, linewidth = 0, label = "data histogram")
histogram!(rand(mix_travel,length(x)), nbins = 100, fillalpha = 0.3, linewidth = 0, xlims = (0,60),label = "simulation histogram")
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/Travel_Mix_Gamma_simulation_Cas.pdf")
##
lims = (0,0.135)
plot(0:60,pdf(mix_travel,0:60),ylims = lims, linecolor = :red, label = "Convex combination")
plot!(0:60,pdf(mix_travel.components[2],0:60),ylims = lims, linecolor = :cyan, xticks = 0:5:60, label = "First component")
plot!(0:60,pdf(mix_travel.components[1],0:60),ylims = lims, linecolor = :magenta, label = "Second component")
q95 = quantile(mix_travel.components[2],0.95)
vline!([q95], label = "95th percentile\nfirst component")
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/Travel_Mix_Gamma_selection_criteria_Cas.pdf")
##
cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    r.AfterLast <= 5 &&
    r.Travel_to < q95
    ,pre_cas)
cas_travel = dvAnalysis(cas_df,:Virus,:Travel_to)
cas_travel.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/TravelTimeCas.pdf")
##
q05 = quantile(mix_travel.components[1],0.05)
cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    r.AfterLast <= 5 &&
    r.Travel_to > q05
    ,pre_cas)
cas_travel = dvAnalysis(cas_df,:Virus,:Travel_to)
cas_travel.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/SecondTravelTimeCas.pdf")
## Travel time Age
pre_age = filter(r -> 0 < r.Travel_to <= 60, Age_s)
pre_age[!,:Travel_to] = round.(pre_age.Travel_to, digits = 2)
x = pre_age.Travel_to
p = @df pre_age histogram(:Travel_to, nbins = 100,label = "data histogram")
m = p.series_list[1].plotattributes[:y]

mix_travel = mixture_gamma(x)
scaling = maximum(filter(!isnan,m))/maximum(pdf(mix_travel,1:60))
plot!(0:60,pdf(mix_travel,0:60)*scaling, label = "mixture distribution fit", linecolor = :red)
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/Travel_Mix_Gamma_fit_Juv.pdf")
##
@df pre_age histogram(:Travel_to, nbins = 100, color = :grey, fillalpha = 0.3, linewidth = 0, label = "data histogram")
histogram!(rand(mix_travel,length(x)), nbins = 100, fillalpha = 0.3, linewidth = 0, label = "simulation histogram",xlims = (-1,60))
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/Travel_Mix_Gamma_simulation_Juv.pdf")
##
lims = (0,0.25)
plot(0:60,pdf(mix_travel,0:60),ylims = lims, linecolor = :red, label = "Convex combination")
plot!(0:60,pdf(mix_travel.components[2],0:60),ylims = lims, linecolor = :cyan, xticks = 0:5:60, label = "First component")
plot!(0:60,pdf(mix_travel.components[1],0:60),ylims = lims, linecolor = :magenta, label = "Second component")
q95 = quantile(mix_travel.components[2],0.95)
vline!([q95], label = "95th percentile\nfirst component")
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/Travel_Mix_Gamma_selection_criteria_Juv.pdf")
##
ag_df = filter(r->
    r.AfterLast <= 5 &&
    r.Travel_to < q95
    ,pre_age)
age_travel = dvAnalysis(ag_df,:Age,:Travel_to)
age_travel.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/TravelTimeJuv.pdf")
##
q05 = quantile(mix_travel.components[1],0.05)
ag_df = filter(r->
    r.AfterLast <= 5 &&
    r.Travel_to > q05
    ,pre_age)
age_travel = dvAnalysis(ag_df,:Age,:Travel_to)
age_travel.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/SecondTravelTimeJuv.pdf")
##
df = agedf
grouping = [:Age,:Sex]
vary = :AfterLast
df1 = individual_summary(df, grouping, vary)
df2 = group_summary(df1, grouping, vary)
df2 = combine(groupby(df1,grouping)) do dd
    m = mean(dd[:,vary])
    s = sem(dd[:,vary])
    (Central = m, ERR = (s,s))
end

FLPDevelopment.mean_sem_scatter(df,grouping,vary)[2]
FLPDevelopment.mode_ci_scatter(df,grouping,vary)
println(unique(agedf[:,:MouseID]))
#########################
quantile(Age_s.AfterLast,0.95)
Age_s.NL = Age_s.AfterLast .<= 7
FLPDevelopment.frequency_plot(Age_s, :NL, :AfterLast)
gd = groupby(Age_s,[:MouseID,:Limit])
df1 = combine(gd, :AfterLast => mean=> :AfterLast)
open_html_table(df1)
df2 = combine(groupby(df1,:Limit)) do dd
    m = mean(dd[:,:AfterLast])
    s = sem(dd[:,:AfterLast])
    (Central = m, ERR = (s,s))
end
FLPDevelopment.frequency_plot(Cas_s, :Limit, :AfterLast)
gd = groupby(Cas_s,[:MouseID,:Limit])
df1 = combine(gd, :AfterLast => maximum => :Max_AfterLast, :AfterLast => length => :Trials)
df1 = combine(gd, :AfterLast => (x-> (union(x),)) => :Max_AfterLast, :AfterLast => length => :Trials)
open_html_table(df1)
df2 = combine(groupby(df1,:Limit)) do dd
    m = mean(dd[:,:AfterLast])
    s = sem(dd[:,:AfterLast])
    (Central = m, ERR = (s,s))
end
###############
combine(groupby(agedf, :Age), :MouseID => x -> length(union(x)))
combine(groupby(casdf, :Virus), :MouseID => x -> length(union(x)))
##
kcas = group_kde(Cas_s,:Trial_duration, group = :Virus)
@df kcas plot(:Xaxis,:Mean, ribbon = :Sem, group = :Virus)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","Caspase.pdf"))

kage = group_kde(Age_s,:Trial_duration, group = :Age)
@df kage plot(:Xaxis,:Mean, ribbon = :Sem, group = :Age, xlims = (0,300))
