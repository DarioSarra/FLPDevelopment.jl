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
## Projection selection
cas_df = filter(r->
    r.Gen == "Rbp4-cre" &&
    r.ProtocolSession == 1  &&
    r.AfterLast <= quantile(Cas_s.AfterLast,0.90),Cas_s)
gdc = groupby(cas_df,[:Virus,:MouseID])
c1 = combine(:AfterLast => mean => :AfterLast,gdc)
c1[!,:xpos] = [v == "tdTomato" ? 1 : 2  for v in c1.Virus]
crecas = filter(r -> r.Virus == "Caspase",c1).AfterLast
cretd = filter(r -> r.Virus == "tdTomato",c1).AfterLast
## test normality to decide between parametric or non parametric tests
ExactOneSampleKSTest(crecas,Normal(mean(crecas),std(crecas)))
ExactOneSampleKSTest(cretd,Normal(mean(cretd),std(cretd)))
## calculate parametric stats
c2 = combine(groupby(c1,:Virus)) do df
    m = mean(df.AfterLast)
    ci = (m - confint(OneSampleTTest(df.AfterLast))[1], confint(OneSampleTTest(df.AfterLast))[2] - m)
    (Mean = mean(df.AfterLast), ERR = ci)
end
@df c2 scatter(1:nrow(c2),:Mean, yerror = :ERR,
    xlims = (0.5, nrow(c2) + 0.5),
    xticks = (1:nrow(c2),:Virus),
    legend = false)
@df c1 scatter!(:xpos,:AfterLast, markersize = 3, alpha = 0.5, color = :grey)
pvalue(UnequalVarianceTTest(crecas,cretd))
## calculate nonparametric stats)
c3 = combine(groupby(c1,:Virus)) do df
    wci = confint(SignedRankTest(collect(df.AfterLast)))
    med = median(df.AfterLast)
    (Median = med, WERR = (med - wci[1],wci[2] - med))
end
@df c3 scatter(1:nrow(c3),:Median, yerror = :WERR,
    xlims = (0.5, nrow(c3) + 0.5),
    xticks = (1:nrow(c3),:Virus),
    legend = false)
@df c1 scatter!(:xpos,:AfterLast, markersize = 3, alpha = 0.5, color = :grey)
pvalue(MannWhitneyUTest(crecas,cretd))
## Age selection
age_df = filter(r->
    r.ProtocolSession == 1 &&
    r.AfterLast <= quantile(Cas_s.AfterLast,0.90),
    Age_s)
gda = groupby(age_df,[:Gen,:MouseID])
a1 = combine(:AfterLast => mean => :AfterLast,gda)
sort!(a1,:Gen)
a1[!,:xpos] = [g == "Adult" ? 1 : 2  for g in a1.Gen]
young = filter(r -> r.Gen == "Young",a1).AfterLast
adult = filter(r -> r.Gen == "Adult",a1).AfterLast
## test normality to decide between parametric or non parametric tests
ExactOneSampleKSTest(young,Normal(mean(young),std(young)))
ExactOneSampleKSTest(cretd,Normal(mean(cretd),std(cretd)))
## calculate parametric stats
a2 = combine(groupby(a1,:Gen)) do df
    m = mean(df.AfterLast)
    ci = (m - confint(OneSampleTTest(df.AfterLast))[1], confint(OneSampleTTest(df.AfterLast))[2] - m)
    (Mean = mean(df.AfterLast), ERR = ci)
end
sort!(a2,:Gen)
@df a2 scatter(1:nrow(a2),:Mean, yerror = :ERR,
    xlims = (0.5, nrow(a2) + 0.5),
    xticks = (1:nrow(a2),:Gen),
    legend = false)
@df a1 scatter!(:xpos,:AfterLast, markersize = 3, alpha = 0.5, color = :grey)
pvalue(UnequalVarianceTTest(crecas,cretd))
## calculate nonparametric stats
a3 = combine(groupby(a1,:Gen)) do df
    wci = confint(SignedRankTest(collect(df.AfterLast)))
    med = median(df.AfterLast)
    (Median = med, WERR = (med - wci[1],wci[2] - med))
end
@df a3 scatter(1:nrow(c3),:Median, yerror = :WERR,
    xlims = (0.5, nrow(c3) + 0.5),
    # xticks = nothing,
    xticks = (1:nrow(c3),:Gen),
    legend = false)
@df a1 scatter!(:xpos,:AfterLast, markersize = 3, alpha = 0.5, color = :grey)
pvalue(MannWhitneyUTest(crecas,cretd))
##
dvplot(cas_df,:Virus,:AfterLast)
dvplot(age_df,:Gen,:AfterLast)
