using Revise, FLPDevelopment, BrowseTables

Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Test"
Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Dev_raw_data"
##
DataIndex, p, b, s = process_dataset(Directory)
DataIndex = get_DataIndex(Directory)
######### Filtering Datas #############
data = filter(r->r.ExpSession == 1,s)
data[!,:Age] = [x in youngs ? "Young" : "Adult" for x in data.MouseID]
data[!,:Exp_type] = [x in age_exp ? "Development" : "Projections" for x in data.MouseID]
data[!,:Treatment] = [x in cnos_animals ? "CNO" : "SAL" for x in data.MouseID]
data[!,:Gen] = [x in wt_animals ? "WT" : "HET" for x in data.MouseID]
data[!,:Combo] = [g*t for (g,t) in zip(data.Gen, data.Treatment)]
data[!,:LowPerformance] = [x in below_thrs for x in data.MouseID]
# filter!(r-> !(r.MouseID in below_thrs), data)
CSV.write("/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Dev_raw_data/Results_2020-05-24/full_streaks.csv",data)
proj = filter(r -> r.Exp_type == "Projections",data)
CSV.write("/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Dev_raw_data/Results_2020-05-24/cno_streaks.csv",proj)
######### Analysis Projection Experiment #########
#data filtering
replace!(proj.PreInterpoke, missing => 0.0) # missing to 0 nothing to worry this data will excluded anyway beacuse it occurs only if they poke just one time
filt = filter(r -> r.Num_pokes > 1 && # at least 2 pokes per trial
    r.Streak > 10 && #exclude the first 10 trials when they have no idea what is going on
    r.PreInterpoke <= 10 && # exclude streaks where the interpoke is longer than 10 seconds
    !r.LowPerformance, # exclude animlas that did less than 15 trials in their session
    proj);
### this is for plotting the mean plus sem bar plot ###
df1 = combine(:AfterLast => mean, groupby(filt,[:Combo, :MouseID]))
df2 = combine(:AfterLast_mean => a -> (Mean = mean(a), SEM = sem(a)),groupby(df1,:Combo))
@df df2 bar(:Combo,:Mean, yerror  = :SEM,
    yticks = 0:0.25:4,
    linecolor = :black, markerstrokecolor = :black)
savefig("/home/beatriz/Documents/Dario_reports/development/mean_sem.pdf");
## this is to plot the boxplot ##
@df filt boxplot(:Combo,:AfterLast)
savefig("/home/beatriz/Documents/Dario_reports/development/boxplot.pdf");
## here is the analysis explained in the email
hetcno = filter(r -> r.Combo == "HETCNO",filt).AfterLast;
hetsal = filter(r -> r.Combo == "HETSAL",filt).AfterLast;
wtcno = filter(r -> r.Combo == "WTCNO",filt).AfterLast;
pvalue(KSampleADTest(hetcno,hetsal,wtcno))
pvalue(MannWhitneyUTest(hetcno,hetsal))
pvalue(MannWhitneyUTest(hetcno,wtcno))
## Alternative analysis
pvalue(KruskalWallisTest(wtcno,hetsal,wtcno))
pvalue(KruskalWallisTest(hetcno,hetsal))
pvalue(KruskalWallisTest(wtcno,hetsal))
