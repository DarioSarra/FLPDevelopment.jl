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
    filter!(r -> r.Protocol == "90/90" &&
    r.ProtocolSession == 1
    ,df)
end
##
#=
 Young vs Adults mice are in the Age labelled dataset (Age_p: pokes, Age_b: bouts Age_s: streaks)
 Caspase vs tdTomato mice are in the Cas labelled dataset (Cas_p: pokes, Cas_b: bouts Cas_s: streaks)
=#
## Afterlast df selection
cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    # r.P_AfterLast >= 0.06
    r.AfterLast <= 5
    # r.AfterLast_frequency >= 0.05
    ,Cas_s)
age_df = filter(r->
    # r.P_AfterLast >= 0.06
    r.AfterLast <= 5
    # r.AfterLast_frequency >= 0.05
    ,Age_s)
maximum(age_df.AfterLast)
## Afterlast Plots
cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast; yspan = (1,2.5))
# cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast)
cas_afterlast.plot
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/AfterLastCas.pdf")
age_afterlast = dvAnalysis(age_df,:Age,:AfterLast; yspan = (1,2.5))
# age_afterlast = dvAnalysis(age_df,:Age,:AfterLast)
age_afterlast.plot
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/AfterLastJuv.pdf")
## Probability df selection
cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    r.ProtocolSession == 1
    ,Cas_s)

age_df = filter(r->
    r.ProtocolSession == 1
    ,Age_s)
## Probability Plots
cas_correct = dvAnalysis(cas_df,:Virus,:CorrectLeave,yspan = (0,1))
cas_correct.plot
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/CorrectCas.pdf")
age_correct = dvAnalysis(age_df,:Age,:CorrectLeave, yspan = (0,1))
age_correct.plot
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/CorrectJuv.pdf")

## Pokes  example
cases = findall((Cas_p.Reward .== false) .& (Cas_p.Correct .== true) .& (Cas_p.PokeInStreak .== 2))
##
idx = cases[4]
p = plot(;legend = false)
# from row 13 to 19
for i in idx - 3:idx + 3
    poke_plot!(p,Cas_p[i,:])
end
p

## Probability mass function Caspase
cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    r.ProtocolSession == 1
    ,Cas_s)
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
df4 = filter(r -> r.AfterLast <= limit, df3)
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
## Travel time Caspase
pre_cas = filter(r -> 1 <= r.Travel_to <= 60, Cas_s)
pre_cas[!,:Travel_to] = round.(pre_cas.Travel_to, digits = 2)
mix_travel = mixture_gamma(pre_cas.Travel_to)
p = @df pre_cas histogram(:Travel_to, nbins = 100,label = "data histogram")
m = p.series_list[1].plotattributes[:y]
scaling = maximum(filter(!isnan,m))/maximum(pdf(mix_travel,0:60))
plot!(0:60,pdf(mix_travel,0:60)*scaling, label = "mixture distribution fit")
##
@df pre_cas histogram(:Travel_to, nbins = 100, color = :grey, fillalpha = 0.3, linewidth = 0, label = "data histogram")
histogram!(rand(mix_travel,length(x)), nbins = 100, fillalpha = 0.3, linewidth = 0, label = "simulation histogram")
##
lims = (0,0.135)
plot(0:60,pdf(mix_travel,0:60),ylims = lims, linecolor = :magenta, label = "Convex combination")
plot!(0:60,pdf(mix_travel.components[2],0:60),ylims = lims, linecolor = :cyan, xticks = 0:5:60, label = "First component")
plot!(0:60,pdf(mix_travel.components[1],0:60),ylims = lims, linecolor = :red, label = "Second component")
q95 = quantile(mix_travel.components[2],0.95)
vline!([q95], label = "95th percentile\nfirst component")
##
cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    r.AfterLast <= 5 &&
    2 <= r.Travel_to < q95
    ,pre_cas)
cas_travel = dvAnalysis(cas_df,:Virus,:Travel_to)
cas_travel.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/TravelTimeCas.pdf")
##
## Travel time Age
pre_age = filter(r -> 1 <= r.Travel_to <= 60, Age_s)
pre_age[!,:Travel_to] = round.(pre_age.Travel_to, digits = 2)
mix_travel = mixture_gamma(pre_age.Travel_to)

p = @df pre_age histogram(:Travel_to, nbins = 100,label = "data histogram")
m = p.series_list[1].plotattributes[:y]
scaling = maximum(filter(!isnan,m))/maximum(pdf(mix_travel,0:60))

plot!(0:60,pdf(mix_travel,0:60)*scaling, label = "mixture distribution fit")
plot(0:60,pdf(mix_travel,0:60), label = "mixture distribution fit")

##
@df pre_age histogram(:Travel_to, nbins = 100, color = :grey, fillalpha = 0.3, linewidth = 0, label = "data histogram")
histogram!(rand(mix_travel,length(x)), nbins = 100, fillalpha = 0.3, linewidth = 0, label = "simulation histogram")
##
lims = (0,0.135)
plot(0:60,pdf(mix_travel,0:60),ylims = lims, linecolor = :magenta, label = "Convex combination")
plot!(0:60,pdf(mix_travel.components[2],0:60),ylims = lims, linecolor = :cyan, xticks = 0:5:60, label = "First component")
plot!(0:60,pdf(mix_travel.components[1],0:60),ylims = lims, linecolor = :red, label = "Second component")
q95 = quantile(mix_travel.components[2],0.95)
vline!([q95], label = "95th percentile\nfirst component")
##
cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    r.AfterLast <= 5 &&
    2 <= r.Travel_to < q95
    ,pre_age)
cas_travel = dvAnalysis(cas_df,:Virus,:Travel_to)
cas_travel.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/TravelTimeCas.pdf")
pre_age = filter(r -> 1 <= r.Travel_to <= 60, Age_s)
pre_age[!,:Travel_to] = round.(pre_age.Travel_to, digits = 2)
@df pre_age histogram(:Travel_to, bins = 100, xticks = 0:5:60, grid = true)
limit = quantile(pre_age.Travel_to,0.95)
age_df = filter(r->
    r.AfterLast <= 5 &&
    r.Travel_to < 7
    ,pre_age)
age_travel = dvAnalysis(age_df,:Age,:Travel_to, yspan = (0,9))
age_travel.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/TravelTimeJuv.pdf")
## Interpoke interval time
limit = quantile(collect(skipmissing(Cas_p.PreInterpoke)),0.95)
cas_df = filter(r ->
    r.Gen == "HET" &&
    !ismissing(r.PreInterpoke) &&
    r.PreInterpoke < limit
    ,Cas_p)
cas_interpoke = dvAnalysis(cas_df,:Virus,:PreInterpoke; yspan = (0,6))
cas_interpoke.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/InterpokeCas.pdf")
##
limit = quantile(collect(skipmissing(Age_p.PreInterpoke)),0.95)
age_df = filter(r ->
    !ismissing(r.PreInterpoke) &&
    r.PreInterpoke < limit
    ,Age_p)
age_interpoke = dvAnalysis(age_df,:Age,:PreInterpoke; yspan = (0,6))
age_interpoke.plot
##
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/InterpokeJuv.pdf")
