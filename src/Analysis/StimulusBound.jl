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
include("Filters.jl")
##
contrasts = Dict(
    :Num_Rewards => StandardizedPredictors.Center(0),
    :Age => DummyCoding(; base="Adults"),
    :Virus => DummyCoding(; base = "tdTomato"),
    :MouseID => Grouping())
## Test if interaction has more explanatory power
f_add_SB_Age = @formula(AfterLast ~ 1 + Num_Rewards+Age +
    (1|MouseID) + (Num_Rewards|MouseID));
f_mult_SB_Age = @formula(AfterLast ~ 1 + Num_Rewards*Age +
    (1|MouseID) + (Num_Rewards|MouseID));
add_SB_Age = fit(MixedModel,f_add_SB_Age,  Age_s; contrasts)
mult_SB_Age = fit(MixedModel,f_mult_SB_Age,  Age_s; contrasts)
MixedModels.likelihoodratiotest(add_SB_Age,mult_SB_Age)


f_add_SB_Cas = @formula(AfterLast ~ 1 + Num_Rewards+Virus +
    (1|MouseID) + (Num_Rewards|MouseID));
f_mult_SB_Cas = @formula(AfterLast ~ 1 + Num_Rewards*Virus +
    (1|MouseID) + (Num_Rewards|MouseID));
add_SB_Cas = fit(MixedModel,f_add_SB_Cas,  Cas_s; contrasts)
mult_SB_Cas = fit(MixedModel,f_mult_SB_Cas,  Cas_s; contrasts)
MixedModels.likelihoodratiotest(add_SB_Cas,mult_SB_Cas)
## plot gentle slope
plot_xy(Age_s,:Num_Rewards,:AfterLast; group = :Age, bin = false, legend = :top,
    ylims = (0,15), ylabel = "Consecutive failures", xlabel = "Rewards obtained")

plot_xy(Cas_s,:Num_Rewards,:AfterLast; group = :Virus, bin = false, legend = :top,
    ylims = (0,15), ylabel = "Consecutive failures", xlabel = "Rewards obtained")
## plot individual coefficients
m = mult_SB_Cas
main_ef = table_coef(m)
ran_ef = table_ranef(m)
ran_ef[!,:Virus] = [get(VirusDict, m,"NA") for m in ran_ef.MouseID]
ineraction_eff = main_ef[4,2]
ran_ef[:, Symbol("Net_Num_Rewards(centered: 0)")] = [v == "tdTomato" ? c : c + ineraction_eff
    for (v,c) in zip(ran_ef.Virus, ran_ef[:,Symbol("Net_Num_Rewards(centered: 0)")])]
@df ran_ef scatter(:Virus, cols(Symbol("Net_Num_Rewards(centered: 0)")))
##
stimbound_f = @formula(AfterLast ~ 1 + Num_Rewards  + Age + (1|MouseID)+ (Num_Rewards|MouseID));
stimbound_m = fit(MixedModel,stimbound_f, Age_s; contrasts)
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
cas_stimbound_f = @formula(AfterLast ~ 1 + Num_Rewards  + Virus + (1|MouseID)+ (Num_Rewards|MouseID));
cas_stimbound_m = fit(MixedModel,cas_stimbound_f, Cas_s; contrasts)

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
Age_SB = fit(MixedModel,f_Age_SB,  Age_s; contrasts)
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
Cas_SB = fit(MixedModel,f_Cas_SB,  Cas_s; contrasts)
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
