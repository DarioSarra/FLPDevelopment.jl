using Revise, FLPDevelopment, BrowseTables
gr(size=(600,600), tick_orientation = :out, grid = false,
    linecolor = :black,
    markerstrokecolor = :black,
    thickness_scaling = 2,
    markersize = 8)
##
include("Young_to_run.jl")
include("Caspase_to_run.jl")
for df in (Age_p, Age_b, Age_s, Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Protocol == "90/90", df)
end
Cas_s[!,:Gen] = [g == "HET" ? "Rbp4-cre" : "Wild Type" for g in Cas_s.Gen]
Cas_s[!,:Combo] = Cas_s.Gen .* "\n" .* Cas_s.Virus
Cas_s[!,:Group] = [r.Gen == "Wild Type" ? "Wild Type" : r.Gen .* "\n" .* r.Virus for r in eachrow(Cas_s)]
categorical!(Cas_s,[:Gen,:Virus,:Group])
levels!(Cas_s.Gen,["Wild Type","Rbp4-cre",])
levels!(Cas_s.Virus,["tdTomato", "Caspase"])
levels!(Cas_s.Group,["Rbp4-cre\ntdTomato", "Rbp4-cre\nCaspase", "Wild Type",])
##
#=
 Young vs Adults mice are in the Age labelled dataset (Age_p: pokes, Age_b: bouts Age_s: streaks)
 Caspase vs tdTomato mice are in the Cas labelled dataset (Cas_p: pokes, Cas_b: bouts Cas_s: streaks)
=#
##
cas_df = filter(r->
    r.Gen == "Rbp4-cre" &&
    r.ProtocolSession == 1  &&
    r.AfterLast <= quantile(Cas_s.AfterLast,0.90),Cas_s)
gd1 = groupby(cas_df,[:Virus,:MouseID])
c1 = combine(:AfterLast => mean => :AfterLast,gd1)
c1[!,:xpos] = [v == "tdTomato" ? 1 : 2  for v in c1.Virus]
c2 = combine(groupby(c1,:Virus)) do df
    (Mean = mean(df.AfterLast), ERR = confint(OneSampleTTest(df.AfterLast))[2] - mean(df.AfterLast))
end
@df c2 scatter(1:nrow(c2),:Mean, yerror = :ERR,
    xlims = (0.5, nrow(c2) + 0.5),
    xticks = (1:nrow(c2),:Virus),
    legend = false)
@df c1 scatter!(:xpos,:AfterLast, markersize = 3, alpha = 0.5, color = :grey)
crecas = filter(r -> r.Virus == "Caspase",c1).AfterLast
cretd = filter(r -> r.Virus == "tdTomato",c1).AfterLast
pvalue(EqualVarianceTTest(crecas,cretd))
pvalue(MannWhitneyUTest(crecas,cretd))
pvalue(ApproximateTwoSampleKSTest(crecas,cretd))
##
age_df = filter(r->
    r.ProtocolSession == 1 &&
    r.AfterLast <= quantile(Cas_s.AfterLast,0.90),
    Age_s)
gd1 = groupby(age_df,[:Gen,:MouseID])
c1 = combine(:AfterLast => mean => :AfterLast,gd1)
c1[!,:xpos] = [g == "Adult" ? 2 : 1  for g in c1.Gen]
c2 = combine(groupby(c1,[:Gen])) do df
    (Mean = mean(df.AfterLast), ERR = confint(OneSampleTTest(df.AfterLast))[2] - mean(df.AfterLast))
end

@df c2 scatter(1:nrow(c2),:Mean, yerror = :ERR,
    xlims = (0.5, nrow(c2) + 0.5),
    xticks = (1:nrow(c2),:Gen),
    legend = false)
@df c1 scatter!(:xpos,:AfterLast, markersize = 3, alpha = 0.5, color = :grey)
young = filter(r -> r.Gen == "Young",c1).AfterLast
adult = filter(r -> r.Gen == "Adult",c1).AfterLast
pvalue(EqualVarianceTTest(young,adult))
pvalue(MannWhitneyUTest(young,adult))
pvalue(ApproximateTwoSampleKSTest(young,adult))
