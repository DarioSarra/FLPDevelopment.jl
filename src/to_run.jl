using Revise, FLPDevelopment, BrowseTables

Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Test"
Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Dev_raw_data"
DataIndex, p, b, s = process_dataset(Directory)
######### Filtering Datas #############
data = filter(r->r.ExpSession == 1,s)
data[!,:Age] = [x in youngs ? "Young" : "Adult" for x in data.MouseID]
data[!,:Exp_type] = [x in age_exp ? "Development" : "Projections" for x in data.MouseID]
data[!,:Treatment] = [x in cnos_animals ? "CNO" : "SAL" for x in data.MouseID]
data[!,:Gen] = [x in wt_animals ? "WT" : "HET" for x in data.MouseID]
data[!,:Combo] = [g*t for (g,t) in zip(data.Gen, data.Treatment)]
gd = groupby(data,:MouseID)
transform!(gd, :Streak => maximum => :Performance)
CSV.write("/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Dev_raw_data/Results_2020-05-24/full_streaks.csv",data)
proj = filter(r -> r.Exp_type == "Projections",data)
CSV.write("/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Dev_raw_data/Results_2020-05-24/cno_streaks.csv",proj)
######### Analysis Projection Experiment #########
#data filtering
@df proj scatter(:MouseID, :Performance, group = :Gen,xrotation = 50, yticks = 0:10:300,
    xlabel = "Mouse ID", ylabel = "Trials number")
savefig("/home/beatriz/Documents/Dario_reports/development/performance.pdf");
replace!(proj.PreInterpoke, missing => 0.0) # missing to 0 nothing to worry this data will excluded anyway beacuse it occurs only if they poke just one time
filt = filter(r -> r.Num_pokes > 1 && # at least 2 pokes per trial
    r.Streak > 10 && #exclude the first 10 trials when they have no idea what is going on
    r.PreInterpoke <= 10 && # exclude streaks where the interpoke is longer than 10 seconds
    r.Performance > 49, # exclude animlas that did less than 15 trials in their session
    proj);
### this is for plotting the mean plus sem bar plot ###
df1 = combine(:AfterLast => mean, groupby(filt,[:Combo, :MouseID]))
df2 = combine(:AfterLast_mean => a -> (Mean = mean(a), SEM = sem(a)),groupby(df1,:Combo))
@df df2 bar(:Combo,:Mean, yerror  = :SEM,
    yticks = 0:0.25:4,
    linecolor = :black, markerstrokecolor = :black,
    xlabel = " Gen x Treatment", ylabel = "Pokes After Last Rew, mean + sem")
savefig("/home/beatriz/Documents/Dario_reports/development/50Trials_mean_sem.pdf");
## this is to plot the boxplot ##
@df filt boxplot(:Combo,:AfterLast,xlabel = " Gen x Treatment", ylabel = "Pokes After Last Rew, median + CI")
savefig("/home/beatriz/Documents/Dario_reports/development/50Trials_boxplot.pdf");
## Pokes after last reward stats per trial
hetcno = filter(r -> r.Combo == "HETCNO",filt).AfterLast;
hetsal = filter(r -> r.Combo == "HETSAL",filt).AfterLast;
wtcno = filter(r -> r.Combo == "WTCNO",filt).AfterLast;
pvalue(KSampleADTest(hetcno,hetsal,wtcno))
pvalue(MannWhitneyUTest(hetcno,hetsal))
pvalue(MannWhitneyUTest(wtcno,hetsal))
## Alternative analysis
pvalue(KruskalWallisTest(wtcno,hetsal,wtcno))
pvalue(KruskalWallisTest(hetcno,hetsal))
pvalue(KruskalWallisTest(wtcno,hetsal))
##
df1 = combine([:AfterLast,:Performance, :Gen] => (a,p,g) ->
    (AfterLast = mean(a), Performance = p[1], Gen = g[1]),
    groupby(filt,[:Combo, :MouseID]))
@df df1 boxplot(:Combo,:AfterLast,
    xlabel = " Gen x Treatment",
    ylabel = "Mean After Last Rew by Mouse, median + CI")
savefig("/home/beatriz/Documents/Dario_reports/development/50Trials_boxplotbymouse.pdf");
## Mean pokes after last reward by mouse stats
hetcno = filter(r -> r.Combo == "HETCNO",df1).AfterLast;
hetsal = filter(r -> r.Combo == "HETSAL",df1).AfterLast;
wtcno = filter(r -> r.Combo == "WTCNO",df1).AfterLast;
pvalue(KSampleADTest(hetcno,hetsal,wtcno))
pvalue(MannWhitneyUTest(hetcno,hetsal))
pvalue(MannWhitneyUTest(hetcno,wtcno))
##
@df df1 scatter(:MouseID, :Performance, group = :Combo)
savefig("/home/beatriz/Documents/Dario_reports/development/50Trials_performances.pdf");
################# fraction of errors analysis ###############
gd = groupby(filt, [:MouseID,:Combo])
df1 = combine(:CorrectLeave => x -> (FractionErrors = (1 - length(findall(x))/length(x)) * 100,),gd)
@df df1 boxplot(:Combo,:FractionErrors,
    xlabel = " Gen x Treatment",
    ylabel = "Incorrect leaving percentage by mouse, median + CI" )
savefig("/home/beatriz/Documents/Dario_reports/development/50Trials_Boxplot_errors%.pdf");
##
df2 = combine(:FractionErrors => a -> (Mean = mean(a), SEM = sem(a)),groupby(df1,:Combo))
@df df2 bar(:Combo,:Mean, yerror  = :SEM,
    yticks = 0:0.25:4,
    linecolor = :black, markerstrokecolor = :black,
    xlabel = " Gen x Treatment",
    guidefont = 8,
    ylabel = "Incorrect leaving percentage by mouse \n mean + sem" )
savefig("/home/beatriz/Documents/Dario_reports/development/50Trials_Mean_errors%.pdf");
## Analysis of fraction of errors per trial by mouse
hetcno = filter(r -> r.Combo == "HETCNO",df1).FractionErrors;
hetsal = filter(r -> r.Combo == "HETSAL",df1).FractionErrors;
wtcno = filter(r -> r.Combo == "WTCNO",df1).FractionErrors;
pvalue(KSampleADTest(hetcno,hetsal,wtcno))
pvalue(MannWhitneyUTest(hetcno,hetsal))
pvalue(MannWhitneyUTest(wtcno,hetsal))
