using Revise, FLPDevelopment, BrowseTables
gr(size=(600,600), tick_orientation = :out, grid = false,
    linecolor = :black,
    markerstrokecolor = :black,
    thickness_scaling = 1,
    markersize = 8);
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
## Densities
## No filters
filt_Age_p = filter(r -> 0 < r.PreInterpoke, Age_p)
a = check_distributions(Age_s,filt_Age_p)
maintitle!(a[1], "Age No filters")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age1_No_filtersViolin.pdf"))
maintitle!(a[2], "Age No filters")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age2_No_filtersViolin.pdf"))
filt_Cas_p = filter(r -> 0 < r.PreInterpoke, Cas_p)
c = check_distributions(Cas_s, filt_Cas_p)
maintitle!(c[1],"Virus No filters")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas1_No_filtersViolin.pdf"))
maintitle!(c[2],"Virus No filters")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas2_No_filtersViolin.pdf"))
ALage = group_frequency(Age_s,:AfterLast; group = :Age)
@df ALage plot(:Xaxis, :Mean, yerror = :Sem, group = :Age, linecolor = :auto)
## 30s Duration
f_Age_s,f_Age_p = FLPDevelopment.joinfilter(Age_s,Age_p,:Trial_duration, 30)
filt_Age_p = filter(r -> 0 < r.PreInterpoke, filt_Age_p)
maintitle!(check_distributions(f_Age_s,f_Age_p),"Age 30s Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age_30sDuration.pdf"))
f_Cas_s,f_Cas_p = FLPDevelopment.joinfilter(Cas_s,Cas_p,:Trial_duration, 30)
filt_Cas_p = filter(r -> 0 < r.PreInterpoke, f_Cas_p)
maintitle!(check_distributions(f_Cas_s,f_Cas_p),"Virus 30s Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas_30sDuration.pdf"))
## q95AfterLast
f_Age_s,f_Age_p = FLPDevelopment.joinfilter(Age_s,Age_p,:AfterLast, 7)
filt_Age_p = filter(r -> 0 < r.PreInterpoke, filt_Age_p)
maintitle!(check_distributions(f_Age_s,f_Age_p),"Age q95 Afterlast")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age_q95AfterLast.pdf"))
f_Cas_s,f_Cas_p = FLPDevelopment.joinfilter(Cas_s,Cas_p,:AfterLast, 7)
filt_Cas_p = filter(r -> 0 < r.PreInterpoke, f_Cas_p)
maintitle!(check_distributions(f_Cas_s,f_Cas_p),"Virus q95 Afterlast")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas_q95AfterLast.pdf"))
## 30s Duration q95AfterLast
f_Age_s,f_Age_p = FLPDevelopment.joinfilter(Age_s,Age_p,:Trial_duration, 30)
f_Age_s,f_Age_p = FLPDevelopment.joinfilter(f_Age_s,f_Age_p,:AfterLast, 7)
maintitle!(check_distributions(f_Age_s,f_Age_p),"Age 30s Duration q95 Afterlast")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age_q95AfterLast30sDuration.pdf"))
f_Cas_s,f_Cas_p = FLPDevelopment.joinfilter(Cas_s,Cas_p,:Trial_duration, 30)
f_Cas_s,f_Cas_p = FLPDevelopment.joinfilter(Cas_s,Cas_p,:AfterLast, 7)
maintitle!(check_distributions(f_Cas_s,f_Cas_p),"Virus 30s Duration q95 Afterlast")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas_q95AfterLast30sDuration.pdf"))
## PreInterpoke
filt_Age_s, filt_Age_p = FLPDevelopment.process_filtered_streak(Age_p,:PreInterpoke,1)
filt_Age_p = filter(r -> 0 < r.PreInterpoke, filt_Age_p)
maintitle!(check_distributions(filt_Age_s,filt_Age_p),"Age Interpoke < 1")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age_Interpoke<1.pdf"))
filt_Cas_s, filt_Cas_p =FLPDevelopment.process_filtered_streak(Cas_p,:PreInterpoke,1)
filter!(r -> 0 < r.PreInterpoke, filt_Cas_p)
maintitle!(check_distributions(filt_Cas_s,filt_Cas_p),"Virus Interpoke < 1")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas_Interpoke<1.pdf"))
## No filters using median for times
filt_Age_p = filter(r -> 0 < r.PreInterpoke, Age_p)
maintitle!(check_distributions(Age_s,filt_Age_p; summary_opt = :MEDIAN), "Age No filters with median")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age_No_filters_Median.pdf"))
filt_Cas_p = filter(r -> 0 < r.PreInterpoke, Cas_p)
maintitle!(check_distributions(Cas_s, filt_Cas_p; summary_opt = :MEDIAN),"Virus No filters with median")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas_No_filters_Median.pdf"))
############# Analysis
ALAge0 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak  + (1|MouseID)),Age_s))
ALAge1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Age + (1|MouseID)),Age_s))
Likelyhood_Ratio_test(ALAge0,ALAge1)

ALAge00 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),Age_s,Poisson())
ALAge01 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Age + (1|MouseID)),Age_s,Poisson())
Likelyhood_Ratio_test(ALAge00,ALAge01)

ALCas0 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak  + (1|MouseID)),Cas_s))
ALCas1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),Cas_s))
Likelyhood_Ratio_test(ALCas0,ALCas1)

ALCas00 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),Cas_s,Poisson())
ALCas01 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),Cas_s,Poisson())
Likelyhood_Ratio_test(ALCas00,ALCas01)
########################## Using each mouse 95% quantile of AfterLast as filter #######################
agedf = filter(r -> r.Limit, Age_s)
casdf = filter(r -> r.Limit, Cas_s)

ALAge0 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak  + (1|MouseID)),agedf))
ALAge1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Age + (1|MouseID)),agedf))
ALAge2 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Age + Sex + (1|MouseID)),agedf))
Likelyhood_Ratio_test(ALAge0,ALAge1)
Likelyhood_Ratio_test(ALAge1,ALAge2)
DoubleAnalysis(agedf,:Age,:AfterLast)

ALAge00 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),agedf,Poisson())
ALAge01 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Age + (1|MouseID)),agedf,Poisson())
Likelyhood_Ratio_test(ALAge00,ALAge01)

ALCas0 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak  + (1|MouseID)),casdf))
ALCas1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),casdf))
Likelyhood_Ratio_test(ALCas0,ALCas1)

ALCas00 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),casdf,Poisson())
ALCas01 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),casdf,Poisson())
Likelyhood_Ratio_test(ALCas00,ALCas01)

filt_Cas_p = filter(r -> 0 < r.PreInterpoke, Cas_p)
c = check_distributions(casdf, filt_Cas_p)
maintitle!(c[1],"Virus No filters")

p_df = group_frequency(casdf,:AfterLast; group = :Virus)
filter!(r -> !isnan(r.Sem), p_df)
@df p_df plot(:Xaxis, :Mean, ribbon = :Sem, linecolor = :auto, group = :Virus)

FLPDevelopment.vars_fig(casdf)
df1 = individual_summary(casdf, :Virus, :AfterLast)

##
casdf = filter(r -> r.Streak <= 100, Cas_s)
check = DoubleAnalysis(Cas_s,:Virus,:AfterLast)
check.nonparametric_plot
check2 = DoubleAnalysis(casdf,:Virus,:AfterLast)
check2.nonparametric_plot
##
#AfterLast Density per mouse
p = plot()
for m in union(Cas_s.MouseID)
    ongoing = filter(r -> r.MouseID == m, Cas_s)
    p_df = FLPDevelopment.single_frequency(ongoing,:AfterLast)
    @df p_df plot!(:Xaxis, :Vals, label = m, linecolor = :auto)
end
p
open_html_table(Age_s)
gd = groupby(Age_s,[:Session,:MouseID])
combine(gd) do dd
    dd1 = filter(r -> !r.Limit, dd)
    [union(dd1.AfterLast)]
end
##
mdf = DataFrame!(CSV.File("/home/beatriz/Desktop/odour_dose_data_FM - Sheet1.csv"))
mdf[!, :Stim] = categorical([x == 0 for x in mdf.Stim])
categorical(mdf.Stim)
hdf = filter(r -> r.Genotype =="sert", mdf)
base = fit(MixedModel,@formula(BOLD ~ 1 + Concentration + (1|MouseID)),hdf)
alt1 = fit(MixedModel,@formula(BOLD ~ 1 + Concentration + Stim + (1|MouseID)),hdf)
Likelyhood_Ratio_test(base,alt1)
AICc_test(base, alt1)
alt2 = fit(MixedModel,@formula(BOLD ~ 1 + Concentration * Stim + (1|MouseID)),hdf)
Likelyhood_Ratio_test(alt1,alt2)
alt3 = fit(MixedModel,@formula(BOLD ~ 1 + Concentration&Stim + (1|MouseID)),hdf)
AICc_test(base, alt3)
