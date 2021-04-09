using Revise, FLPDevelopment, BrowseTables, Random
plotlyjs(size=(600,600), tick_orientation = :out, grid = false,
    # linecolor = :black,
    markerstrokecolor = :black,
    thickness_scaling = 1,
    markersize = 8)
##
include("Young_to_run2.jl")
include("Caspase_to_run.jl")
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
for df in (Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Gen == "Rbp4-cre", df)
end
open_html_table(FLPDevelopment.summarydf(Age_s,Age_p))
open_html_table(FLPDevelopment.summarydf(Cas_s,Cas_p))
##
##
# Logistic regression of the leave choice Leave ~ Streak + PokeIn + PokeHierarchy + Group
# zscore full dataset all variable
# logscale poke In before model
# with predict or with the data
# bin In
# by mouse mean P of leave per bin
# mean group + sem per animal
############### Leaving Model
# zscore full dataset all variable
# logscale poke Out before model?
Cas_p.LogOut = log10.(Cas_p.Out)
transform!(Cas_p, [:Streak, :Out, :LogOut] .=> zscore)
transform!(Cas_s, :Streak .=> zscore)

#Compare AfterLast model
verbagg00 = @formula(Leave ~ 1 + Streak_zscore + Out_zscore + (1|MouseID));
LeaveCas00 = fit(MixedModel,verbagg00, Cas_p, Bernoulli())
verbagg0 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + (1|MouseID));
LeaveCas0 = fit(MixedModel,verbagg0, Cas_p, Bernoulli())
AlCas0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak_zscore + (1|MouseID)),Cas_s)
AIC_test(AlCas0, LeaveCas0)
AIC_test(LeaveCas00, LeaveCas0)

verbagg1 = @formula(Leave ~ 1 + Streak_zscore * Virus + LogOut_zscore * Virus +  (1|MouseID));
verbagg2 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore * Virus +  (1|MouseID));
verbagg3 = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + Virus + (1|MouseID));

# Group Effect
LeaveCas1 = fit(MixedModel,verbagg1, Cas_p, Bernoulli())
LeaveCas2 = fit(MixedModel,verbagg2, Cas_p, Bernoulli())
LeaveCas3 = fit(MixedModel,verbagg3, Cas_p, Bernoulli())


MixedModels.likelihoodratiotest(LeaveCas0,LeaveCas1)
MixedModels.likelihoodratiotest(LeaveCas0,LeaveCas3)


Cas_p.M1_Leave = predict(LeaveCas1)
Cas_p.M2_Leave = predict(LeaveCas2)
Cas_p.M3_Leave = predict(LeaveCas3)
#### Row data P leaving
transform!(Cas_p, :Out => (x -> bin_axis(x; unit_step = 10)) => :Bin_Out)
df0 = filter(r -> r.Streak <= 50, Cas_p)
vars = [:Leave, :M1_Leave, :M2_Leave]
df1 = combine(groupby(df0,[:MouseID,:Bin_Out,:Virus]),
    vars .=> mean .=> vars)
df2 = combine(groupby(df1,[:Bin_Out,:Virus]),
    vars .=> mean .=> [:Mean_leave, :Mean_M1Leave, :Mean_M2Leave],
    vars .=> sem .=> [:Sem_leave, :Sem_M1Leave, :Sem_M2Leave])
filter!(r -> !isnan(r.Sem_leave), df2)
sort!(df2,:Bin_Out)
@df df2 plot(:Bin_Out, :Mean_leave, ribbon = :Sem_leave, group = :Virus,
    linecolor = :auto, xlims = (5,80), ylims = (0,1), legend = :top,
    xlabel = "Poke time from trial beginning (s)",
    ylabel = "Probability of leaving")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","1-Cas_Pleave.pdf"))

transform!(Cas_p, :LogOut => (x -> bin_axis(x; length = 20)) => :Bin_LogOut)
df0 = filter(r -> r.Streak <= 50, Cas_p)
df1 = combine(groupby(df0,[:MouseID,:Bin_LogOut,:Virus]),
    vars .=> mean .=> vars)
df2 = combine(groupby(df1,[:Bin_LogOut,:Virus]),
    vars .=> mean .=> [:Mean_leave, :Mean_M1Leave, :Mean_M2Leave],
    vars .=> sem .=> [:Sem_leave, :Sem_M1Leave, :Sem_M2Leave])
filter!(r -> !isnan(r.Sem_leave), df2)
sort!(df2,:Bin_LogOut)
@df df2 plot(:Bin_LogOut, :Mean_leave, ribbon = :Sem_leave, group = :Virus,
    linecolor = :auto, ylims = (0,1),
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Probability of leaving")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","2-Cas_LogPleave.pdf"))
## PokeN over trials
Cas_s.Num_pokes
Cas_s[!,:BinnedStreak] = bin_axis(Cas_s.Streak; unit_step = 5)
res = summary_xy(Cas_s,:BinnedStreak,:AfterLast; group = :Virus)
cas_alXtrial = @df res plot(string.(:BinnedStreak),:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, xlabel = "Trial", ylabel = "Pokes after last reward")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","2-Cas_AlTrial.pdf"))
res2 = summary_xy(Cas_s,:BinnedStreak,:Num_pokes; group = :Virus)
cas_pokesXtrial = @df res2 plot(string.(:BinnedStreak),:Mean, group = :Virus, linecolor = :auto,
    ribbon = :Sem, xrotation = 50, xlabel = "Trial", ylabel = "Number of pokes",size=(600,600))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","2-Cas_pokesTrial.pdf"))
## Model Plot
rng = MersenneTwister(1234321)
samp1 = parametricbootstrap(rng,100,LeaveCas1)
df = DataFrame(samp1.allpars)
bootdf = combine(groupby(df,[:type, :group, :names]), :value => shortestcovint => :interval)
bootdf.coef = push!(coef(LeaveCas1), mean(ranef(LeaveCas1)[1]))
bootdf.variable = ["Intercept", "Trial", "Virus:Caspase", "Poke-time",
    "Trial & Virus: Caspase",
    "Poke-time & Virus:Caspase", "MouseID"]
transform!(bootdf, [:coef, :interval] => ByRow((c,e) -> (c -e[1], e[2]-c)) => :err)
@df bootdf[1:end-1,:] scatter(:coef ,1:nrow(bootdf)-1,
    xerror = :err, xlabel = "Coefficient estimate",
    yticks = (1:nrow(bootdf), :variable))
vline!([0], linecolor = :red, legend = false)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","3-LogModel.pdf"))
## P leave Heatmaps
df1 = copy(Cas_p)
transform!(df1, :LogOut => (x -> bin_axis(x; length = 30)) => :Bin_LogOut)#7
transform!(df1, :Streak => (x -> bin_axis(x; unit_step = 10)) => :Bin_Streak)
vars = [:Leave, :M1_Leave, :M2_Leave, :M3_Leave]
df2 = combine(groupby(df1,[:Virus,:Bin_LogOut,:Bin_Streak]),
    vars .=> mean .=> vars)
##
LeaveType = :Leave
dfcas = filter(r -> r.Virus == "Caspase", df2)
sort!(dfcas,[:Bin_Streak,:Bin_LogOut])
filter!(r -> -0.4 <= r.Bin_LogOut <= 2.1 && r.Bin_Streak <= 71, dfcas)
reshapecas = unstack(dfcas,  :Bin_Streak, :Bin_LogOut, LeaveType)
ylab = string.(sort(union(reshapecas.Bin_Streak)))
xlab = string.(sort(union(dfcas.Bin_LogOut)))
heatmap(xlab, ylab, Matrix(reshapecas[:, Not(:Bin_Streak)]),
    clim = (0, 1),
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
##savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","4-HeatmapCasRaw.pdf"))
# tom
dftom = filter(r -> r.Virus == "tdTomato", df2)
sort!(dftom,[:Bin_Streak,:Bin_LogOut])
filter!(r -> -0.4 <= r.Bin_LogOut <= 2.1 && r.Bin_Streak <= 71, dftom)
reshapetom = unstack(dftom,:Bin_Streak, :Bin_LogOut, LeaveType)
ylab = string.(sort(union(reshapetom.Bin_Streak)))
xlab = string.(sort(union(dftom.Bin_LogOut)))
heatmap(xlab, ylab, Matrix(reshapetom[:, Not(:Bin_Streak)]),
    clim = (0, 1),
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
##savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","5-HeatmapTomRaw.pdf"))
# difference
differenza = Matrix(reshapecas[:, Not(:Bin_Streak)]) - Matrix(reshapetom[:, Not(:Bin_Streak)])
heatmap(xlab, ylab, differenza,
    clim = (-1, 1),
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
##savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","6-HeatmapDiffRaw.pdf"))
##
minOut = 0.1
maxOut = 1.9
minStreak = 11
maxStreak = 71
LeaveType = :Leave
dfcas = filter(r -> r.Virus == "Caspase", df2)
sort!(dfcas,[:Bin_Streak,:Bin_LogOut])
filter!(r -> minOut <= r.Bin_LogOut <= maxOut && minStreak <= r.Bin_Streak <= maxStreak, dfcas)
reshapecas = unstack(dfcas,  :Bin_Streak, :Bin_LogOut, LeaveType)
ycas = string.(sort(union(reshapecas.Bin_Streak)))
xcas = string.(sort(union(dfcas.Bin_LogOut)))
heatmap(xcas, ycas, Matrix(reshapecas[:, Not(:Bin_Streak)]),
    clim = (0, 1),
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
##savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","4-ZoomHeatmapCasRaw.pdf"))
# tom
dftom = filter(r -> r.Virus == "tdTomato", df2)
sort!(dftom,[:Bin_Streak,:Bin_LogOut])
filter!(r -> minOut <= r.Bin_LogOut <= maxOut &&  minStreak <= r.Bin_Streak <= maxStreak, dftom)
reshapetom = unstack(dftom,:Bin_Streak, :Bin_LogOut, LeaveType)
ytom = string.(sort(union(reshapetom.Bin_Streak)))
xtom = string.(sort(union(dftom.Bin_LogOut)))
heatmap(xtom, ytom, Matrix(reshapetom[:, Not(:Bin_Streak)]),
    clim = (0, 1),
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
##savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","5-ZoomHeatmapTomRaw.pdf"))
# difference
differenza = Matrix(reshapecas[:, Not(:Bin_Streak)]) - Matrix(reshapetom[:, Not(:Bin_Streak)])
ydiff = string.(sort(union(reshapetom.Bin_Streak)))
xdiff = string.(sort(union(dftom.Bin_LogOut)))
heatmap(xdiff, ydiff, differenza,
    clim = (-1, 1),
    xlabel = "Poke time from trial beginning (log10 s)",
    ylabel = "Trial",
    colorbar_title = "P Leave",
    color = :deep)
##savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Leaving","6-ZoomHeatmapDiffRaw.pdf"))
