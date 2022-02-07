using Revise, FLPDevelopment, BrowseTables
## Plot settings
# Font size and type
xyfont = font(18, "Bookman Light")
legfont = font(14, "Bookman Light")

# Figure/Plot image size
pixel_factor = 600
ps_x = pixel_factor * 1
ps_y = pixel_factor * 1
fig_size = (ps_x, ps_y)

theme(:default)
gr(size = fig_size,
    titlefont = font(14, "Bookman Light"),
    guidefont = font(14, "Bookman Light"),
    tick_orientation = :out,
    grid = false,
    markerstrokecolor = :black,
    markersize = 8,
    thickness_scaling = 1.5,
    legendfont = font(10, "Bookman Light"))
# pgfplotsx(size = fig_size,
#     tick_orientation = :out,
#     grid = false,
#     markerstrokecolor = :black,
#     markersize = 8,
#     thickness_scaling = 1.5,
#     titlefont = xyfont,
#     guidefont = xyfont,
#     legendfont = legfont)
## Load and adjust data
include("Young_to_run2.jl")
include("Caspase_to_run.jl")
for df in (Age_p, Age_b, Age_s, Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Protocol == "90/90" &&
    r.MouseID != "CD09" && # biting, see B1_CD09_2020-07-13 minute30
    r.MouseID != "RJ58" && # blind
    r.MouseID != "RJ67" && # biting, see B3_RJ67_2020-09-28 minute 7:33
    !(r.MouseID in first_females_group) &&
    r.ProtocolSession == 1
    ,df)
end
for df in (Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Gen == "Rbp4-cre", df) # exclude wild type animals
end
#adjustments for pokes DataFrames
for pokedf in [Age_p, Cas_p]
    # transform time in log 10 scale
    pokedf.LogOut = log10.(pokedf.Out)
    pokedf.LogInterpoke = log10.(pokedf.PostInterpoke)
    pokedf.LogDuration = log10.(pokedf.PokeDur)
    # zscore values for logistic regression
    transform!(pokedf, [:Streak, :Out, :LogOut] .=> zscore)
end
#adjustments for trials DataFrames
for streakdf in [Age_s, Cas_s]
    # transform time in log 10 scale
    streakdf.LogDuration = log10.(streakdf.Trial_duration)
    # streakdf.ElapsedTime = log10.(streakdf.AfterLast_Duration .+ 0.1)
    streakdf.LogTravel = log10.(streakdf.Travel_to)
    # bin trials
    # streakdf.BinnedStreak = bin_axis(streakdf.Streak; unit_step = 4)
    # streakdf.BlockBinnedStreak = bin_axis(streakdf.Streak; unit_step = 15)
    streakdf.LongBinnedStreak = bin_axis(streakdf.Streak; unit_step = 20)
    transform!(streakdf, :Streak .=> zscore)
end
open_html_table(FLPDevelopment.summarydf(Age_s,Age_p))
open_html_table(FLPDevelopment.summarydf(Cas_s,Cas_p))
########################### Example Session Plots ######################################
## Adult example session
lowertrial = 20
uppertrial = 45
tt = filter(r -> r.MouseID == "RJ23" && lowertrial <= r.Streak <= uppertrial, Age_p)
tt = filter(r -> r.MouseID == "RJ23", Age_p)
t1 = filter(r -> r.Leave, tt)
sort!(t1,[:Out])
order = Dict(t1.Streak .=> 1:nrow(t1))
tt.order = [get(order,x,0) for x in tt.Streak]
sort!(tt, :order)
plt = session_plot(tt; ordered = true)
xaxis!(plt, xlims = :auto)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Pokes","AdultSession.pdf"))
## Juvenile example session
tt = filter(r -> r.MouseID == "RJ02" && lowertrial <= r.Streak <= uppertrial, Age_p)
tt = filter(r -> r.MouseID == "RJ02", Age_p)
t1 = filter(r -> r.Leave, tt)
t2 = sort(t1,[:Out])
order = Dict(1:nrow(t2) .=> t2.Streak)
tt.order = [get(order,x,0) for x in tt.Streak]
sort!(tt, :order)
plt = session_plot(tt)
xaxis!(plt, xlims = :auto)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Pokes","JuvenileSession.pdf"))
## Well trained example session
# comp_df = CSV.read("/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/Datasets/Stimulations/DRN_Opto_again/pokesDRN_Opto_again.csv", DataFrame)
# comp_df = DataFrame(CSV.File(joinpath("/Volumes/GoogleDrive/My Drive/Flipping/Datasets/","Stimulations/DRN_Opto_again/pokesDRN_Opto_again.csv")))
mother1 = replace(path,basename(path)=>"")[1:end-1]
mother2 = replace(mother1,basename(mother1)=>"")
comp_df = CSV.read(joinpath(mother2,"Stimulations","DRN_Opto_again","pokesDRN_Opto_again.csv"))
comp90 = filter(r-> r.Protocol == "90/90", comp_df)
sessions = union(comp90.Session)
tt = filter(r-> r.Session == sessions[99] && lowertrial <= r.Streak <= uppertrial, comp90)
transform!(groupby(tt, [:Session, :Streak]), [:PokeIn, :PokeOut] => ((i,o) -> (In = i .- i[1], Out = o .- i[1])) => AsTable)
session_plot(tt)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Fig2","TrainedSession.pdf"))
## Density leaving times
Age = Age_s[:,[:MouseID,:Trial_duration]]
Age[!,:Condition] = "Naive_".*Age_s[:,:Age]
comp90 = filter(r-> r.Protocol == "90/90", comp_df)
comp_df = CSV.read(joinpath(mother2,"Stimulations","DRN_Opto_again","streaksDRN_Opto_again.csv"))
Comp = comp90[:,[:MouseID,:Trial_duration]]
Comp[!,:Condition] .=  "Trained"
whole = vcat(Age,Comp)
whole[!,:LogDuration] = log10.(whole.Trial_duration)
@df whole density(:LogDuration, group = :Condition,
    xlabel = "Elapsed time (log10 seconds)", ylabel = "Kernell density estimate")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Fig2","DensityTrialDuration.pdf"))
## Age filtered trial duration
Age_max60 = filter(r->r.Trial_duration <= 60, Age_s)
Age_filt = Difference(Age_max60,:Age,:Trial_duration, ind_summary = mean,
    ylabel = "Trial duration (seconds)", ylims = (0,19))
Age_filt.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Fig2","Age_FilteredTrialDur.pdf"))
## Age Survival Plot
filter(r-> r.Age == "Juveniles",Age_s)
filter(r-> r.Virus == "tdTomato",Cas_s)

Age_LevingAn, Age_LevingAn_df = function_analysis(Age_s,:LogDuration, cumulative_algorythm; grouping = :Age, calc = :bootstrapping)
xprop = ("Poke Time(seconds)", xyfont,(log10.([0.1,1,10,100,1000]),["0.1","1","10","100","1000"]))
yprop = ("Probablity of leaving", xyfont)
plot!(Age_LevingAn, xaxis = xprop, yaxis = yprop, legend = false)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","LevingRate.pdf"))
##
Age_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Basic = fit(MixedModel,Age_Basic_verb, Age_p, Bernoulli())
Age_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Age + LogOut_zscore * Age +  (1+Streak_zscore+LogOut_zscore|MouseID));
Age_Full = fit(MixedModel,Age_Full_verb, Age_p, Bernoulli())
MixedModels.likelihoodratiotest(Age_Basic,Age_Full)
Age_noTrial_verb = @formula(Leave ~ 1 + LogOut_zscore * Age +  (1+LogOut_zscore|MouseID));
Age_noTrial =  fit(MixedModel,Age_noTrial_verb, Age_p, Bernoulli())
MixedModels.likelihoodratiotest(Age_noTrial,Age_Full)
Age_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000AgeBootstrap.csv"),DataFrame)
# Age_BootDf = bootstrapdf(Age_p, Age_Full; n = 1000)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","1000AgeBootstrap.csv"),Age_BootDf)
Age_btdf = Age_BootDf[[1,4,2,3,5,6],:]
Age_btdf.err = [tuple(parse.(Float64,split(err[2:end-1],", "))...) for err in Age_btdf.err]
Age_btdf.variable
yprop = ("",font(10, "Bookman Light"),(collect(1:nrow(Age_btdf)),
    [L"Intercept",L"PokeTime",L"Trial",L"Juveniles",
    L"Trial \&",
    L"PokeTime \&"]))
    xprop = ("Coefficient estimate", xyfont, (-1.33,1.3))
    Titolo = L"+ \Leftarrow Behavioural Control \Rightarrow -"
    @df Age_btdf scatter(:coef ,1:nrow(Age_btdf), xerror = :err,legend = false,
    xaxis = xprop, yaxis = yprop, markercolor = :gray75, title = Titolo)
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
## Num Pokes
Age_np = Difference(Age_s, :Age, :Num_pokes; ylabel = "Number of pokes", ylims = (0,6.5))
Age_np.plt
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Pokes","AgeNpokes.pdf"))
Cas_np = Difference(Cas_s, :Virus, :Num_pokes; ylabel = "Number of pokes", ylims = (0,6.5))
Cas_np.plt
Cas_np.test
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Pokes","CasNpokes.pdf"))
## Caspase example session
tt = filter(r -> r.MouseID == "CD10", Cas_p)
combine(groupby(Cas_p, [:MouseID]), :Virus => union)
t1 = filter(r -> r.Leave, tt)
sort!(t1,[:Out])
order = Dict(t1.Streak .=> 1:nrow(t1))
tt.order = [get(order,x,0) for x in tt.Streak]
sort!(tt, :order)
plt = session_plot(tt; ordered = true)
xaxis!(plt, xlims = :auto)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Pokes","CaspaseSession.pdf"))
## tdTomato example session
tt = filter(r -> r.MouseID == "CD17", Cas_p)
combine(groupby(Cas_p, [:MouseID]), :Virus => union)
t1 = filter(r -> r.Leave, tt)
sort!(t1,[:Out])
order = Dict(t1.Streak .=> 1:nrow(t1))
tt.order = [get(order,x,0) for x in tt.Streak]
sort!(tt, :order)
plt = session_plot(tt; ordered = true)
xaxis!(plt, xlims = :auto)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Pokes","tdTomatoSession.pdf"))
## Virus trial duration density
@df Cas_s density(:LogDuration, group = :Virus,
    xlabel = "Elapsed time (log10 seconds)", ylabel = "Kernell density estimate")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","December2021","Fig2","VirusDensityTrialDuration.pdf"))
