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
##
#=
 Young vs Adults mice are in the Age labelled dataset (Age_p: pokes, Age_b: bouts Age_s: streaks)
 Caspase vs tdTomato mice are in the Cas labelled dataset (Cas_p: pokes, Cas_b: bouts Cas_s: streaks)
=#
## Afterlast df selection

cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    r.ProtocolSession == 1 &&
    r.P_AfterLast >= 0.06
    # r.AfterLast <= 5
    # r.AfterLast_frequency >= 0.05
    ,Cas_s)
maximum(cas_df.AfterLast)
age_df = filter(r->
    r.ProtocolSession == 1 &&
    r.P_AfterLast >= 0.06
    # r.AfterLast <= 5
    # r.AfterLast_frequency >= 0.05
    ,Age_s)
maximum(age_df.AfterLast)
## Afterlast Plots
cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast; yspan = (0,2.5))
# cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast)
cas_afterlast.plot
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/AfterLastCas.pdf")
age_afterlast = dvAnalysis(age_df,:Age,:AfterLast; yspan = (0,2.5))
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
gdc = groupby(cas_df,[xvar,:MouseID])
df1 = combine(yvar => mean => yvar,gdc)
sort!(df1,xvar)
df2 = combine(groupby(df1,[xvar])) do dd
    m = mean(dd[:,yvar])
    SEM = sem(dd[:,yvar])
    (Central = m, ERR = SEM)
end
df3 = filter(r -> r.AfterLast <= 25, df2)
df3
@df df3 scatter(:AfterLast,:Central,
    # group = :Virus,
    yerror = :ERR,
    yticks = 0:0.05:0.4, xticks = 0:5:60, grid = true, linecolor = :auto, markersize = 3)
sum(df3.Central)
## Probability mass function Age
age_df = filter(r->
    r.ProtocolSession == 1
    ,Age_s)
xvar = :AfterLast
yvar = :AfterLast_frequency
gdc = groupby(age_df,[xvar,:MouseID])
df1 = combine(yvar => mean => yvar,gdc)
sort!(df1,xvar)
df2 = combine(groupby(df1,[xvar])) do dd
    m = mean(dd[:,yvar])
    SEM = sem(dd[:,yvar])
    (Central = m, ERR = SEM)
end
df3 = filter(r -> r.AfterLast <= 25, df2)
df3
@df df3 scatter(:AfterLast,:Central,
    # group = :Virus,
    yerror = :ERR,
    yticks = 0:0.05:0.4, xticks = 0:5:60, grid = true, linecolor = :auto, markersize = 3)
