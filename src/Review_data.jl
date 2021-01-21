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
    r.MouseID != "CD09" && # biting, see B1_CD09_2020-07-13 minute30
    r.MouseID != "RJ58" && # blind
    # r.MouseID != "RJ67" && # biting, see B3_RJ67_2020-09-28 minute 7:33
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
## Densities
## No filters
check_distributions(Age_s,Age_p)
title!("Age No filters")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age_No_filters.pdf"))
check_distributions(Cas_s, Cas_p)
title!("Virus No filters")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas_No_filters.pdf"))
## 30s Duration 1s Intepoke
f_Age_s = filter(r -> r.Trial_duration < 30, Age_s)
f_Age_p = filter(r -> 0 < r.PreInterpoke < 1, Age_p)
check_distributions(f_Age_s,f_age_p)
title!("Age 30s Duration 1s Intepoke")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age_30sDuration1sIntepokes.pdf"))
f_Cas_s = filter(r -> r.Trial_duration < 30 && r.Gen == "Rbp4-cre", Cas_s)
f_Cas_p = filter(r -> 0 < r.PreInterpoke < 1 && r.Gen == "Rbp4-cre", Cas_p)
check_distributions(f_Cas_s,f_Cas_p)
title!("Virus 30s Duration 1s Intepoke")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas_30sDuration1sIntepokes.pdf"))
## q95AfterLast 1s Intepoke
f_Age_s = filter(r -> r.AfterLast < 7, Age_s)
f_Age_p = filter(r -> 0 < r.PreInterpoke < 1, Age_p)
check_distributions(f_Age_s,f_age_p)
title!("Age q95 Afterlast 1s Intepoke")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age_q95AfterLast1sIntepokes.pdf"))
f_Cas_s = filter(r -> r.AfterLast < 8 && r.Gen == "Rbp4-cre", Cas_s)
f_Cas_p = filter(r -> 0 < r.PreInterpoke < 1 && r.Gen == "Rbp4-cre", Cas_p)
check_distributions(f_Cas_s,f_Cas_p)
title!("Virus q95 Afterlast 1s Intepoke")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas_q95AfterLast1sIntepokes.pdf"))
## 30s Duration q95AfterLast 1s Intepoke
f_Age_s = filter(r -> r.AfterLast < 7 && r.Trial_duration < 30, Age_s)
f_Age_p = filter(r -> 0 < r.PreInterpoke < 1, Age_p)
check_distributions(f_Age_s,f_age_p)
title!("Age 30s Duration q95 Afterlast")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Age_q95AfterLast30sDuration.pdf"))
f_Cas_s = filter(r -> r.AfterLast < 8 && r.Trial_duration < 30 && r.Gen == "Rbp4-cre", Cas_s)
f_Cas_p = filter(r -> 0 < r.PreInterpoke < 1 && r.Gen == "Rbp4-cre", Cas_p)
check_distributions(f_Cas_s,f_Cas_p)
title!("Virus q95 Afterlast 1s Intepoke")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Distributions","Cas_q95AfterLast30sDuration.pdf"))
