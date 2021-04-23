using Revise, FLPDevelopment, BrowseTables
gr(size=(600,600), tick_orientation = :out, grid = false,
    linecolor = :black,
    markerstrokecolor = :black,
    thickness_scaling = 2,
    markersize = 8)
include("Young_to_run2.jl")
include("Caspase_to_run.jl")
##
for df in (Age_p, Age_b, Age_s, Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Protocol == "90/90" &&
    r.MouseID != "CD09" && # biting, see B1_CD09_2020-07-13 minute30
    r.MouseID != "RJ58" && # blind
    r.MouseID != "RJ67" && # biting, see B3_RJ67_2020-09-28 minute 7:33
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
age_al_lim = quantile(Age_s.AfterLast,0.95)
age_df = filter(r->
    # r.Streak <=80
    # r.Performance >25 &&
    r.Trial_duration < 30 &&
    r.AfterLast <= age_al_lim &&
    r.Stop/60 < 50
    ,Age_s)
cas_al_lim = quantile(Cas_s.AfterLast,0.95)
cas_df = filter(r->
    # r.Streak <=80 &&
    # r.Performance >25 &&
    r.Trial_duration < 30 &&
    r.AfterLast < cas_al_lim &&
    r.Stop/60 <= 50  &&
    r.Gen == "Rbp4-cre"
    ,Cas_s)
##
limit = quantile(age_df.Travel_to,0.95)
# limit = maximum(age_df.Travel_to)
shorter = filter(r -> r.Travel_to <= limit, age_df)
x = shorter.Travel_to
mix_dist = mixture_gamma_exp(x)
interval = 0:0.5:80
plot(interval,pdf(mix_dist,interval),xlims = (0,80), label = "Convex combination",
    yrotation = 60)
plt = twinx()
histogram!(plt,x, nbins = 40, color = :grey, fillalpha = 0.3, xlims = (0,80), linewidth = 0, label = false,
    xlabel = "Trial duration")

plot(interval,pdf(mix_dist,interval),xlims = (0,80), label = "Convex combination")
plot(interval,pdf(mix_dist.components[1],interval), xlims = (0,80), linecolor = :cyan, label = "First component")
plot!(interval,pdf(mix_dist.components[2],interval), xlims = (0,80),linecolor = :magenta, label = "Second component")
vline!([quantile(mix_dist.components[2],0.95)])
##
limit = quantile(cas_df.Travel_to,0.95)
# limit = maximum(cas_df.Travel_to)
shorter = filter(r -> r.Travel_to < limit, cas_df)
x = shorter.Travel_to
mix_dist = mixture_gamma(x)
interval = 0:0.5:80
plot(interval,pdf(mix_dist,interval),xlims = (0,80), label = "Convex combination",
    yrotation = 60)
plt = twinx()
histogram!(plt,x, nbins = 40, color = :grey, fillalpha = 0.3, xlims = (0,80), linewidth = 0, label = false,
    xlabel = "Trial duration")

plot(interval,pdf(mix_dist,interval),xlims = (0,80), label = "Convex combination")
plot(interval,pdf(mix_dist.components[1],interval), xlims = (0,80), linecolor = :cyan, label = "First component")
plot!(interval,pdf(mix_dist.components[2],interval), xlims = (0,80),linecolor = :magenta, label = "Second component")
vline!([quantile(mix_dist.components[2],0.95)])
