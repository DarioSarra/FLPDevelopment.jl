using Revise, FLPDevelopment, BrowseTables
gr(size=(600,600), tick_orientation = :out, grid = false,
    linecolor = :black,
    markerstrokecolor = :black,
    thickness_scaling = 2,
    markersize = 8)
##
include("Young_to_run.jl")
include("Caspase_to_run.jl")
##
for df in (Age_p, Age_b, Age_s, Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Protocol == "90/90", df)
end
##
#=
 Young vs Adults mice are in the Age labelled dataset (Age_p: pokes, Age_b: bouts Age_s: streaks)
 Caspase vs tdTomato mice are in the Cas labelled dataset (Cas_p: pokes, Cas_b: bouts Cas_s: streaks)
=#
## Projection selection
cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    r.ProtocolSession == 1 &&
    r.P_AfterLast >= 0.1
    # r.AfterLast <= quantile(Cas_s.AfterLast,0.90)
    ,Cas_s)

cas_afterlast = dvAnalysis(cas_df,:Virus,:AfterLast)
cas_afterlast.plot

cas_afterlast.mice_summary
cas_afterlast.group_summary
cas_afterlast.normality


cas_df = filter(r->
    r.Gen == "Rbp4-cre"&&
    r.ProtocolSession == 1
    # r.AfterLast <= quantile(Cas_s.AfterLast,0.90)
    ,Cas_s)
cas_correct = dvAnalysis(cas_df,:Virus,:CorrectLeave)
cas_correct.plot
## Age selection
age_df = filter(r->
    r.ProtocolSession == 1 &&
    r.P_AfterLast >= 0.1
    # r.AfterLast <= quantile(Cas_s.AfterLast,0.99)
    ,Age_s)
age_afterlast = dvAnalysis(age_df,:Gen,:AfterLast)
age_afterlast.plot



age_df = filter(r->
    r.ProtocolSession == 1
    # r.AfterLast <= quantile(Cas_s.AfterLast,0.99)
    ,Age_s)
age_correct = dvAnalysis(age_df,:Gen,:CorrectLeave)
age_correct.plot
