using Revise, FLPDevelopment, BrowseTables
# Plot settings
# Font size and type
xyfont = font(18, "Bookman Light")
legfont = font(14, "Bookman Light")

# Figure/Plot image size
pixel_factor = 600
ps_x = pixel_factor * 1
ps_y = pixel_factor * 1
fig_size = (ps_x, ps_y)

theme(:default)
# gr(size = fig_size,
#     titlefont = font(14, "Bookman Light"),
#     guidefont = font(14, "Bookman Light"),
#     tick_orientation = :out,
#     grid = false,
#     markerstrokecolor = :black,
#     markersize = 8,
#     thickness_scaling = 1.5,
#     legendfont = font(10, "Bookman Light"))
pgfplotsx(size = fig_size,
    tick_orientation = :out,
    grid = false,
    markerstrokecolor = :black,
    markersize = 8,
    thickness_scaling = 1.5,
    titlefont = xyfont,
    guidefont = xyfont,
    legendfont = legfont)
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
    transform!(groupby(pokedf,[:MouseID, :Streak]), :PokeDur => cumsum => :CumDuration)
    pokedf.Occupancy = pokedf.CumDuration ./ pokedf.Out
    pokedf.Occupancy_zscore = zscore(pokedf.Occupancy)
end
#adjustments for trials DataFrames
for streakdf in [Age_s, Cas_s]
    # transform time in log 10 scale
    streakdf.LogDuration = log10.(streakdf.Trial_Duration)
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
##
stimbound_f = @formula(AfterLast ~ 1 + Num_Rewards  + Age + (1+Num_Rewards|MouseID));
stimbound_m = fit(MixedModel,stimbound_f, Age_s)
simple_age_coeff = DataFrame(only(raneftables(stimbound_m)))
rename!(simple_age_coeff, Symbol("(Intercept)") => :Intercept, :Num_Rewards => :Res_NumRewards)
simple_age_coeff[!,:Coef_NumRewards] = simple_age_coeff.Res_NumRewards .+ stimbound_m.β[2]
simple_age_coeff[!,:Age] = [x in dario_youngs ? "Juveniles" : "Adults" for x in simple_age_coeff.MouseID]
transform!(groupby(simple_age_coeff,:Age), :Intercept =>
    (x-> round.(accumulate(+, repeat([0.02],length(x)); init = 0), digits =2)) => :Shift)
transform!(simple_age_coeff, [:Shift,:Age] =>
    ((p,a) -> round.(p .+ [x == "Juveniles" ? 1 : 2 for x in a], digits = 2)) => :Pos)
plot_xy(Age_s,:Num_Rewards,:AfterLast; group = :Age, bin = false, legend = :top,
    ylims = (0,15), ylabel = "Consecutive failures", xlabel = "Rewards obtained")
extrema(simple_age_coeff.Pos)
@df simple_age_coeff scatter!(:Pos, :Coef_NumRewards, group = :Age,
    markersize = 2, xticks = ([1.25, 2.25], ["Juveniles", "Adults"]), xlims = (1.02 - 0.1 ,2.43 + 0.1),
    legend = false, ylabel = "Correlation coefficients",guidefontsize = 9,
    inset = (1, bbox(0.15, 0.15, 0.4, 0.4, :top, :left)),
    subplot =2)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Post-Rewievs", "StimulusBound.png"))
##
cas_stimbound_f = @formula(AfterLast ~ 1 + Num_Rewards  + Virus + (1+Num_Rewards|MouseID));
cas_stimbound_m = fit(MixedModel,cas_stimbound_f, Cas_s)
simple_cas_coeff = DataFrame(only(raneftables(cas_stimbound_m)))
rename!(simple_cas_coeff, Symbol("(Intercept)") => :Intercept, :Num_Rewards => :Res_NumRewards)
simple_cas_coeff[!,:Coef_NumRewards] = simple_cas_coeff.Res_NumRewards .+ cas_stimbound_m.β[2]
simple_cas_coeff[!,:Virus] = [get(VirusDict,x,"Missing") for x in simple_cas_coeff.MouseID]
transform!(groupby(simple_cas_coeff,:Virus), :Intercept =>
    (x-> round.(accumulate(+, repeat([0.02],length(x)); init = 0), digits =2)) => :Shift)
transform!(simple_cas_coeff, [:Shift,:Virus] =>
    ((p,a) -> round.(p .+ [x == "Caspase" ? 1 : 2 for x in a], digits = 2)) => :Pos)
transform!(simple_cas_coeff, :Virus => categorical => :Virus)
levels!(simple_cas_coeff.Virus,["tdTomato", "Caspase"])
sort!(simple_cas_coeff,[:Virus,:Coef_NumRewards], rev = true)
##
plot_xy(Cas_s,:Num_Rewards,:AfterLast; group = :Virus, bin = false, legend = :topright,
    ylims = (0,6.8), ylabel = "Consecutive failures", xlabel = "Rewards obtained")
extrema(simple_cas_coeff.Pos)
@df simple_cas_coeff scatter!(:Pos, :Coef_NumRewards, group = :Virus,
    markersize = 2, xticks = ([1.1, 2.1], ["Caspase", "tdTomato"]), xlims = (1.02 - 0.1 ,2.3 + 0.1),
    legend = false, ylabel = "Correlation coefficients",guidefontsize = 9,
    inset = (1, bbox(0.1, 0.02, 0.39, 0.39, :top, :left)),
    subplot =2)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Post-Rewievs", "Cas_StimulusBound.png"))
##
Cas_Rewards = Difference(Cas_s, :Virus, :Num_Rewards, ylabel = "number of rewards per trial", ylims = (0,1.2))
Cas_Rewards.plt
Cas_Rewards.test
Cas_Rewards.groupdf
Cas_Rewards.groupdf[!,:Measure] .= "Number of rewards"
rename!(Cas_Rewards.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Post-Rewievs", "Cas_Rewards.pdf"))
##
f_SB = @formula(AfterLast ~ 1 + Num_Rewards + (1+Num_Rewards|MouseID));
Ad_SB = fit(MixedModel,f_SB, filter(r-> r.Age == "Adults", Age_s))
Juv_SB = fit(MixedModel,f_SB, filter(r-> r.Age == "Juveniles", Age_s))
Td_SB = fit(MixedModel,f_SB, filter(r-> r.Virus == "tdTomato", Cas_s))
Cas_SB = fit(MixedModel,f_SB, filter(r-> r.Virus == "Caspase", Cas_s))
##
f_Age_SB = @formula(AfterLast ~ 1 + Num_Rewards*Age + (1+Num_Rewards|MouseID));
Age_SB = fit(MixedModel,f_Age_SB,  Age_s)
age_coeff = DataFrame(only(raneftables(Age_SB)))
rename!(age_coeff, Symbol("(Intercept)") => :Intercept, :Num_Rewards => :mouse_dev_NRew)
age_coeff[!,:Age] = [x in dario_youngs ? "Juveniles" : "Adults" for x in age_coeff.MouseID]
transform!(age_coeff, :Age => ByRow(x -> x == "Adults" ? 0 : Age_SB.β[4]) => :group_dev_NRew)
age_coeff[!,:Coef_NRew] = Age_SB.β[2] .+ (age_coeff.mouse_dev_NRew .+  age_coeff.group_dev_NRew)
transform!(groupby(age_coeff,:Age), :Intercept =>
    (x-> round.(accumulate(+, repeat([0.02],length(x)); init = 0), digits =2)) => :Shift)
transform!(age_coeff, [:Shift,:Age] =>
    ((s,a) -> round.(s .+ [x == "Juveniles" ? 1 : 2 for x in a], digits = 2)) => :Pos)
@df age_coeff scatter(:Pos, :Coef_NRew, group = :Age,
    markersize = 3, xticks = ([1.25, 2.25], ["Juveniles", "Adults"]),
    xlims = (1.02 - 0.1 ,2.43 + 0.1), ylabel = "AfterLast/NumRewards coeficient", legend =  :top)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Post-Rewievs", "Age_coeff.png"))
##
f_Cas_SB = @formula(AfterLast ~ 1 + Num_Rewards*Virus + (1+Num_Rewards|MouseID));
Cas_SB = fit(MixedModel,f_Cas_SB,  Cas_s)
cas_coeff = DataFrame(only(raneftables(Cas_SB)))
rename!(cas_coeff, Symbol("(Intercept)") => :Intercept, :Num_Rewards => :mouse_dev_NRew)
cas_coeff[!,:Virus] = [get(VirusDict,x,"Missing") for x in cas_coeff.MouseID]
transform!(cas_coeff, :Virus => categorical => :Virus)
levels!(cas_coeff.Virus,["tdTomato", "Caspase"])
transform!(cas_coeff, :Virus => ByRow(x -> x == "tdTomato" ? 0 : Cas_SB.β[4]) => :group_dev_NRew)
cas_coeff[!,:Coef_NRew] = Cas_SB.β[2] .+ (cas_coeff.mouse_dev_NRew .+  cas_coeff.group_dev_NRew)
transform!(groupby(cas_coeff,:Virus), :Intercept =>
    (x-> round.(accumulate(+, repeat([0.02],length(x)); init = 0), digits =2)) => :Shift)
transform!(cas_coeff, [:Shift,:Virus] =>
    ((s,a) -> round.(s .+ [x == "Caspase" ? 1 : 2 for x in a], digits = 2)) => :Pos)
@df cas_coeff scatter(:Pos, :Coef_NRew, group = :Virus,
    markersize = 3, xticks = ([1.1, 2.1], ["Caspase", "tdTomato"]),
    xlims = (1.02 - 0.1 ,2.2 + 0.1), ylabel = "AfterLast/NumRewards coeficient", legend = :top)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Post-Rewievs", "Cas_coeff.png"))
##
open_html_table(sort(simple_age_coeff,:Age))
open_html_table(sort(age_coeff, :Age))
open_html_table(sort(simple_cas_coeff,:Virus))
open_html_table(sort(cas_coeff, :Virus))
