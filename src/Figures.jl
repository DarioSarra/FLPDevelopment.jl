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
    r.ProtocolSession == 1
    ,df)
end
##
#=
to do list
- single animal distributions
- group median and ci plus mean animal afterlast plot
- CD9 interpoke, trial duaration and travel time against Rbp4
- CD9 long trials interpoke, trial duaration and travel time against short trials
- CD9 video
=#
#=
 Young vs Adults mice are in the Age labelled dataset (Age_p: pokes, Age_b: bouts Age_s: streaks)
 Caspase vs tdTomato mice are in the Cas labelled dataset (Cas_p: pokes, Cas_b: bouts Cas_s: streaks)
=#
## Afterlast df selection
cas_df = filter(r->
    r.Gen == "Rbp4-cre" &&
    # r.Performance == 25 &&
    r.MouseID != "CD09"
    # r.PreInterpoke <= 1
    # r.AfterLast <= 6
    # r.AfterLast_frequency >= 0.05
    ,Cas_s)
age_df = filter(r->
    # r.P_AfterLast >= 0.06
    r.Sex != "c"
    # r.Performance >= 25
    # r.PreInterpoke <= 1
    # r.AfterLast_frequency >= 0.05
    ,Age_s)
maximum(age_df.AfterLast)
age_df
open_html_table(cas_df)
mouse_summary(cas_df,:Virus,:AfterLast)
########################### Afterlast Plots ######################################
cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast; nonparametric = true, yspan = (0,6))
# cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast)
cas_afterlast.plot
cas_afterlast.normality
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/AfterLastCas.pdf")
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
