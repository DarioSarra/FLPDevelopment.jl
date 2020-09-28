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
- group median and ci plus mean animal afterlast plot
- single animal distributions
- CD9 interpoke, trial duaration and travel time against Rbp4
- CD9 long trials interpoke, trial duaration and travel time against short trials
- CD9 video
=#
## Afterlast df selection
cas_df = filter(r->
    r.Gen == "Rbp4-cre"
    ,Cas_s)
age_df = filter(r->
    r.Sex != "c"
    ,Age_s)
## AfterLast with outlier
cas_afterlast = DoubleAnalysis(cas_df,:Virus,:AfterLast)
cas_afterlast.JarqueBera
cas_afterlast.parametric_plot
cas_afterlast.nonparametric_plot
# savefig("/Volumes/GoogleDrive/My Drive/Reports for Zach/Development project/AfterLastCas.pdf")
## Outlier analysis minute 30:04 starts biting and breaking the beam
o_df = filter(r -> r.Virus == "Caspase",cas_df)
 ###### After Last #######
check_cd9(o_df,:AfterLast)
 ###### InterPoke #######
check_cd9(o_df,:PreInterpoke)
 ###### Travel to #######
check_cd9(o_df,:Travel_to)
##
cd9_s = filter(r -> r.MouseID == "CD09",cas_df)

@df cd9_s scatter(:AfterLast,:PreInterpoke, markersize = 3)
@df cd9_s scatter(:AfterLast,:Travel_to, markersize = 3)

cd9_p = filter(r -> r.MouseID == "CD09",Cas_p)
@df cd9_p scatter(:PokeInStreak,:PreInterpoke, markersize = 3)
open_html_table(cd9_p)
open_html_table(cd9_s)
##
check_cd9(Cas_p,:PokeInStreak,:PreInterpoke)
