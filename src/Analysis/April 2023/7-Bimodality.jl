include(joinpath(splitpath(@__DIR__)[1:end-1]...,"Filters.jl"))
include("Bimodality_tools.jl")
## Cluster LogDuration data in 2 kmeans
kcluster_df!(Age_s,2,:LogDuration)
categorise_kclusters!(Age_s, :LogDuration, :Assignment, ["Short", "Long"]; new_name = :Strategy)
## calculate proportion of strategies per mouse
Age_k = combine(groupby(Age_s,[:Age,:MouseID])) do dd
    diz = countmap(dd.Strategy)
    (Short = diz["Short"]/nrow(dd), Long = diz["Long"]/nrow(dd), Trials = nrow(dd))
end
#calculate difference, but per se we don't expect the difference to be equal to 0
transform!(Age_k, [:Long,:Short] => ByRow((l,s)-> l-s) => :Diff)
Age_group_fit = AkaikeBimodality(Age_k,:Diff; group = :Age) #add params number
Age_mouse_fit = AkaikeBimodality(Age_s,:LogDuration; group = [:Age,:MouseID])
combine(groupby(Age_mouse_fit,:Age), [:Unimodal_AICc, :Bimodal_AICc] => ((u,b) -> sum(b.<u)) => :Bimodals)
## Fisher Exact test
transform!(Age_s, :Strategy => ByRow(x -> x == "Long") => :Long)
Age_f = combine(groupby(Age_s,:Age)) do dd0
    dd2 = DataFrame()
    for m in unique(dd0.MouseID)
        dd_g = filter(r-> r.MouseID != m, dd0)
        dd_m = filter(r-> r.MouseID == m, dd0)
        group_events = sum(dd_g.Long)
        group_observations = length(dd_g.Long)
        mouse_events = sum(dd_m.Long)
        mouse_observations = length(dd_m.Long)
        test = FisherExactTest(group_events, mouse_events, group_observations, mouse_observations)
        dd1 = DataFrame(MouseID = m,P = pvalue(test), CI = confint(test), 
        Mouse_events = mouse_events, Mouse_observations = mouse_observations,
        Group_events = group_events, Group_observations = group_observations,
        P_Long= mouse_events/mouse_observations)
        isempty(dd2) ? (dd2 = dd1) : append!(dd2,dd1)
    end
    return dd2
end
Age_f[!,:Outliers] = Age_f.P .< 0.05
Age_f_grouped = filter(r->r.Outliers,Age_f)
open_html_table(Age_f)
##
combine(groupby(Age_k,:Age), :Diff => (x -> compareGMMs(Vector(x))) => :P_single)
Age_mouse_fit = combine(groupby(Age_s,[:MouseID,:Age]), :LogDuration => (x -> compareGMMs(Vector(x))) => :P_single)
combine(groupby(Age_mouse_fit,:Age), :P_single => (x -> sum(x.<0.05)) => :Bimodals)
## Age Plots
#adjust colors
control_c, manipulation_c = Plots.theme_palette(:auto)[13], Plots.theme_palette(:auto)[14]
col_df = combine(groupby(Age_s,:Age)) do dd
    m = unique(dd.MouseID)
    l = length(m) + 2
    c = dd.Age[1] == "Adults" ? :ice : :amp
    p = palette(c, l)[2:l-1]
    (MouseID  = m, Col = p)
end
col_dic = Dict(id => col for (id,col) in zip(col_df.MouseID,col_df.Col));
Age_s[!,:col] = [get(col_dic,x, :black) for x in Age_s.MouseID];
# plot individual trials strategy assignment
@df sort(Age_s,:Strategy) dotplot(:Strategy, :LogDuration, legend = false,
    markersize = 2, markerstrokewidth = 0, markeralpha = 0.8, color = :col,
    xlabel = "Trial types", #xticks = ([1,2], ["Long", "Short"]),
    yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Leaving time (s)",
    title = "Age manipulation"
)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_ClusteredTrials.png"))
## Plot difference in strategy use
Age_k[!,:col] = [x == "Adults" ? 13 : 14 for x in  Age_k.Age]
@df sort(Age_k,:Age) dotplot(:Age,:Diff, legend = false,
    markersize = 4, markerstrokewidth = 1, color = :col,
    ylims = (-1,1), ylabel = "Difference in strategy frequencies (long - short)",
    title = "Age manipulation"
)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_ClusteredGroup.png"))
## Caspase
kcluster_df!(Cas_s,2,:LogDuration)
categorise_kclusters!(Cas_s, :LogDuration, :Assignment, ["Short", "Long"]; new_name = :Strategy)
## calculate proportion of strategies per mouse
Cas_k = combine(groupby(Cas_s,[:Virus,:MouseID])) do dd
    diz = countmap(dd.Strategy)
    (Short = diz["Short"]/nrow(dd), Long = diz["Long"]/nrow(dd), Trials = nrow(dd))
end
#calculate difference, but per se we don't expect the difference to be equal to 0
transform!(Cas_k, [:Long,:Short] => ByRow((l,s)-> l-s) => :Diff)
Cas_group_fit = AkaikeBimodality(Cas_k,:Diff; group = :Virus) #add params number
Cas_mouse_fit = AkaikeBimodality(Cas_s,:LogDuration; group = [:Virus,:MouseID])
combine(groupby(Cas_mouse_fit,:Virus), [:Unimodal_AICc, :Bimodal_AICc] => ((u,b) -> sum(b.<u)) => :Bimodals)
## Fisher Exact test
transform!(Cas_s, :Strategy => ByRow(x -> x == "Long") => :Long)
Cas_f = combine(groupby(Cas_s,:Virus)) do dd0
    dd2 = DataFrame()
    for m in unique(dd0.MouseID)
        dd_g = filter(r-> r.MouseID != m, dd0)
        dd_m = filter(r-> r.MouseID == m, dd0)
        group_events = sum(dd_g.Long)
        group_observations = length(dd_g.Long)
        mouse_events = sum(dd_m.Long)
        mouse_observations = length(dd_m.Long)
        test = FisherExactTest(group_events, mouse_events, group_observations, mouse_observations)
        dd1 = DataFrame(MouseID = m, P = pvalue(test), CI = confint(test), 
        Mouse_events = mouse_events, Mouse_observations = mouse_observations,
        Group_events = group_events, Group_observations = group_observations,
        P_Long= mouse_events/mouse_observations)
        isempty(dd2) ? (dd2 = dd1) : append!(dd2,dd1)
    end
    return dd2
end
Cas_f[!,:Outliers] = Cas_f.P .< 0.05
Cas_f_grouped = filter(r->r.Outliers,Cas_f)
combine(groupby(Cas_f,:Virus), :Outliers => sum => :Outliers)
open_html_table(Cas_f)
## Cas Plots
#adjust colors
control_c, manipulation_c = Plots.theme_palette(:auto)[13], Plots.theme_palette(:auto)[14]
Cas_coldf = combine(groupby(Cas_s,:Virus)) do dd
    m = unique(dd.MouseID)
    l = length(m) + 2
    c = dd[1, :Virus] == "tdTomato" ? :ice : :amp
    p = palette(c, l)[2:l-1]
    (MouseID  = m, Col = p)
end
col_dic = Dict(id => col for (id,col) in zip(Cas_coldf.MouseID,Cas_coldf.Col));
Cas_s[!,:col] = [get(col_dic,x, :black) for x in Cas_s.MouseID];
# plot individual trials strategy assignment
@df sort(Cas_s,:Strategy) dotplot(:Strategy, :LogDuration, legend = false,
    markersize = 2, markerstrokewidth = 0, markeralpha = 0.8, color = :col,
    xlabel = "Trial types", #xticks = ([1,2], ["Long", "Short"]),
    yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Leaving time (s)",
    title = "Virus manipulation"
)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_ClusteredTrials.png"))
## Plot difference in strategy use
Cas_k[!,:col] = [x == "tdTomato" ? 13 : 14 for x in  Cas_k.Virus]
@df sort(Cas_k,:Virus) dotplot(:Virus,:Diff, legend = false,
    markersize = 4, markerstrokewidth = 1, color = :col,
    ylims = (-1,1), ylabel = "Difference in strategy frequencies (long - short)",
    title = "Virus manipulation"
)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_ClusteredGroup.png"))
open_html_table(Cas_k)
##
out_mouse = ["RJ18", "CD15"]

Age_filt = filter(r -> !in(r.MouseID,out_mouse), Age_p)
Age_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Age_Basic = fit(MixedModel,Age_Basic_verb, Age_filt, Bernoulli(); contrasts)
Age_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Age + LogOut_zscore * Age +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Age_Full = fit(MixedModel,Age_Full_verb, Age_filt, Bernoulli(); contrasts)
MixedModels.likelihoodratiotest(Age_Basic,Age_Full)

Cas_filt = filter(r -> !in(r.MouseID,out_mouse), Cas_p)
Cas_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Cas_Basic = fit(MixedModel,Cas_Basic_verb, Cas_filt, Bernoulli();contrasts)
Cas_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Virus + LogOut_zscore * Virus +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Cas_Full = fit(MixedModel,Cas_Full_verb, Cas_filt, Bernoulli();contrasts)
MixedModels.likelihoodratiotest(Cas_Basic,Cas_Full)
