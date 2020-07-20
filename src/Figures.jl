using Revise, FLPDevelopment, BrowseTables
gr(size=(600,600), tick_orientation = :out, grid = false,
    linecolor = :black,
    markerstrokecolor = :black,
    thickness_scaling = 2)
include("Young_to_run.jl")
include("Caspase_to_run.jl")
for df in (Age_p, Age_b, Age_s, Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Protocol == "90/90", df)
end
#=
 Young vs Adults mice are in the Age labelled dataset (Age_p: pokes, Age_b: bouts Age_s: streaks)
 Caspase vs tdTomato mice are in the Cas labelled dataset (Cas_p: pokes, Cas_b: bouts Cas_s: streaks)
=#
filter!(r->r.AfterLast <= quantile(Age_s.AfterLast,0.95),Age_s)
AfterLast_young = filter(r -> r.Gen == "Young",Age_s).AfterLast
AfterLast_adult = filter(r -> r.Gen == "Adult",Age_s).AfterLast
ApproximateTwoSampleKSTest(AfterLast_young,AfterLast_adult)
res = summarize(Age_s,:Gen,:AfterLast)
@df res scatter([1,2],:Mean, yerror = :SEM,
    xlims = (0.5,2.5),
    xticks = (1:2,:Xaxis))


filter!(r->r.AfterLast <= quantile(Cas_s.AfterLast,0.95),Cas_s)
AfterLast_caspase = filter(r -> r.Virus == "Caspase",Cas_s).AfterLast
AfterLast_tdtomato = filter(r -> r.Virus == "tdTomato",Cas_s).AfterLast
ApproximateTwoSampleKSTest(AfterLast_caspase,AfterLast_tdtomato)
