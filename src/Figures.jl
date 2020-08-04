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
    r.P_AfterLast >= 0.1
    # r.AfterLast <= quantile(Cas_s.AfterLast,0.90)
    ,Cas_s)

age_df = filter(r->
    r.ProtocolSession == 1 &&
    r.P_AfterLast >= 0.1
    # r.AfterLast <= quantile(Cas_s.AfterLast,0.99)
    ,Age_s)
## Afterlast Plots
cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast; yspan = (0,2))
cas_afterlast.plot
savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/AfterLastCas.pdf")
age_afterlast = dvAnalysis(age_df,:Age,:AfterLast; yspan = (0,2))
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
##
findall(.!(Cas_p.Reward) .& Cas_p.Correct)
p = plot(;legend = false)
# from row 13 to 19
for i in 12:19
    poke_plot!(p,Cas_p[i,:])
end
p
##
p = plot(;legend = false)
for i in 403:405#321:333
    poke_plot!(p,Cas_p[i,:])
end
p
##
findall(.!(cas_df.CorrectStart))
i = 52
idx = cas_df[i,:Streak]
ses = cas_df[i,:Session]
df1 = filter(r -> idx-1 <= r.Streak <= idx +1 && r.Session == ses, cas_df)
df2 = filter(r -> idx-3 <= r.Streak <= idx +1 && r.Session == ses, Cas_p)
p = plot(;legend = false)
for r in eachrow(df2)
    poke_plot!(p,r)
end
p
##
cases = findall((Cas_p.Reward .== false) .& (Cas_p.Correct .== true) .& (Cas_p.PokeInStreak .== 2))
##
idx = cases[4]
p = plot(;legend = false)
# from row 13 to 19
for i in idx - 3:idx + 3
    poke_plot!(p,Cas_p[i,:])
end
p
