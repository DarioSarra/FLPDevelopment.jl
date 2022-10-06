include("Filters.jl")
##
contrasts = Dict(
    :Num_Rewards => StandardizedPredictors.Center(0),
    :AfterLast => StandardizedPredictors.Center(0),
    :Age => DummyCoding(; base = "Adults"),
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
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","StimulusBound", "Age_StimBoung.pdf"))
plot_xy(Cas_s,:Num_Rewards,:AfterLast; group = :Virus, bin = false, legend = :top,
    ylims = (0,15), ylabel = "Consecutive failures", xlabel = "Rewards obtained")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","StimulusBound", "Cas_StimBoung.pdf"))

### plot individual coefficients Age
m = mult_SB_Age
main_ef = table_coef(mult_SB_Age)
ran_ef = table_ranef(mult_SB_Age)
ran_ef[!,:Age] = [x in dario_youngs ? "Juveniles" : "Adults" for x in ran_ef.MouseID]
ran_ef[!,:color] = [x == "Adults" ? 1 : 2 for x in ran_ef.Age]
nr_c = main_ef[2,2]
int_c = main_ef[4,2]
v = []
for r in eachrow(ran_ef)
    res = nr_c + r[Symbol("Num_Rewards(centered: 0)")] + int_c * Float64(r.Age == "Juveniles")
    push!(v,res)
end
ran_ef[:, Symbol("Net_Num_Rewards(centered: 0)")] = v
sort!(ran_ef,:Age)
@df ran_ef scatter(:Age, cols(Symbol("Net_Num_Rewards(centered: 0)")),
    xlims = (0,2), legend = false, color = :color, ylabel = "Number of rewards coefficients ")
    hline!([0,0], linecolor = :black, linestyle = :dash)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","StimulusBound", "Age_StimBound_coef.pdf"))
# plot individual coefficients Cas
m = mult_SB_Cas
main_ef = table_coef(m)
ran_ef = table_ranef(m)
open_html_table(main_ef)
    open_html_table(ran_ef)
ran_ef[!,:Virus] = [get(VirusDict, m,"NA") for m in ran_ef.MouseID]
ran_ef[!,:color] = [x == "tdTomato" ? 1 : 2 for x in ran_ef.Virus]
nr_c = main_ef[2,2]
int_c = main_ef[4,2]
v = []
for r in eachrow(ran_ef)
    res = nr_c + r[Symbol("Num_Rewards(centered: 0)")] + int_c * Float64(r.Virus == "Caspase")
    push!(v,res)
end
ran_ef[:, Symbol("Net_Num_Rewards(centered: 0)")] = v
@df ran_ef scatter(:Virus, cols(Symbol("Net_Num_Rewards(centered: 0)")),
    xlims = (0,2), legend = false, color = :color, ylabel = "Number of rewards coefficients ")
    hline!([0,0], linecolor = :black, linestyle = :dash)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","StimulusBound", "Cas_StimBound_coef.pdf"))
##
