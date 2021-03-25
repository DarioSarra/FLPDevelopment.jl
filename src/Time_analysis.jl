using Revise, FLPDevelopment, BrowseTables
gr(size=(600,600), tick_orientation = :out, grid = false,
    linecolor = :black,
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
#####
# raw pokes plot
# #which filtering criteria is a good example
# Travel, Poke, InterPoke, Trial duration: distribution and single animals medians
# histogram short medium long trials
# duration of first poke
# #fraction of poking time in a trial
####

FLPDevelopment.all_duration_analysis(Age_p, Age_s, :Age)[1]
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","AgeTrial.pdf"))
FLPDevelopment.all_duration_analysis(Age_p, Age_s, :Age)[2]
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","AgePoke.pdf"))

FLPDevelopment.all_duration_analysis(Cas_p, Cas_s, :Virus)[1]
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","CasTrial.pdf"))
FLPDevelopment.all_duration_analysis(Cas_p, Cas_s, :Virus)[2]
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","CasPoke.pdf"))

########################### Example Session Plots ######################################
for m in union(Cas_p.MouseID)
    tt = filter(r -> r.MouseID == m, Cas_p)
    plt = plot(legend = false, xlabel = "Trials",
        yticks = 0:10:100,
        ylabel = "Time (seconds)", ylims = (0,100),size=(800,1200))
    for r in eachrow(tt)
        FLPDevelopment.session_plot!(plt, r)
    end
    title!(plt,m * " " * String(tt[1,:Virus]))
    savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime",m*".pdf"))
end

# example biting CD19 minute 3:25
tt = filter(r -> r.MouseID == "CD19" &&  r.Streak == 2, Cas_p)
plt = plot(legend = false, xlabel = "Trials",
    yaxis = false, yticks = false, ylabel = "Time (seconds)",size=(800,800),
    xlims = (150,170), xticks = 0:1:300)
for r in eachrow(tt)
    FLPDevelopment.session_plot!(plt, r)
end
plt
savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","CD19-biting.pdf"))
open_html_table(tt[:,[:PokeInStreak,:PokeDur,:PreInterpoke]])

##
gd = groupby(Cas_s,[:MouseID,:Limit])
df1 = combine(gd, :AfterLast => maximum => :Max_AfterLast, :AfterLast => length => :Trials)
df1 = combine(gd, :AfterLast => (x-> (union(x),)) => :Max_AfterLast, :AfterLast => length => :Trials)
open_html_table(df1)
df2 = combine(groupby(df1,:Limit)) do dd
    m = mean(dd[:,:AfterLast])
    s = sem(dd[:,:AfterLast])
    (Central = m, ERR = (s,s))
end
##
df = Age_s
variable = :Trial_duration
group = :Age
bins = 4
gd = groupby(df, [:MouseID, group])
df2 = combine(gd, variable => (x -> FLPDevelopment.quantile_vec(x, bins)) => variable)
gd2 = groupby(df2, group)
df3 = combine(gd2) do dd
    df3 = DataFrame(Quantile = Float64[], Count = Int[])
    for (k,i) in countmap(dd[:,variable])
        push!(df3, [k,i])
    end
    df3
end
@df df3 groupedbar(:Quantile, :Count, group = cols(group), xticks = 0:1/(bins):1)
###### Binned duration plots
trial = FLPDevelopment.bin_duration(Cas_s,:Trial_duration,:Virus)
plt = plot(trial[1]..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Trial Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusTrialDurations.pdf"))

check = filter(r -> r.Travel_to > 0, Cas_s)
travel = FLPDevelopment.bin_duration(check,:Travel_to,:Virus)
plt = plot(travel[1]..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Travel Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusTravelDurations.pdf"))

check = filter(r -> r.PreInterpoke > 0, Cas_p)
inter = FLPDevelopment.bin_duration(check,:PreInterpoke,:Virus)
check2 = FLPDevelopment.calculate_bin_duration(check,:PreInterpoke,:Virus)
bn = round.(union(check2.Bin), digits =2)
bn = join(string.(bn), " ")
plt = plot(inter[1]..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "InterPoke Duration\n" * bn)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusInterpokeDurations.pdf"))

check = filter(r -> 0 < r.PreInterpoke <= 0.5, Cas_p)
shortinter = FLPDevelopment.bin_duration(check,:PreInterpoke,:Virus; modality = :STEP, unit_step = 0.05)
plt = plot(shortinter[1]..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Short InterPoke Duration\n" * bn)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusShortInterpokeDurations.pdf"))

check = filter(r -> r.Reward, Cas_p)
rew = FLPDevelopment.bin_duration(check,:PokeDur,:Virus)
plt = plot(rew[1]..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Rewards Poke Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusRewardDurations.pdf"))

check = filter(r -> !r.Reward, Cas_p)
failures = FLPDevelopment.bin_duration(check,:PokeDur,:Virus)
plt = plot(failures[1]..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Failures Poke Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusFailuresDurations.pdf"))

check = filter(r -> 0 < !r.Reward && r.PokeDur <= 1, Cas_p)
shortfailures = FLPDevelopment.bin_duration(check,:PokeDur,:Virus; modality = :STEP, unit_step = 0.05)
plt = plot(shortfailures[1]..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Short Failures Poke Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusShortFailuresDurations.pdf"))

plot(
    (@df shortfailures[2] bar(:Bin,:Count_mean, yerr = :Count_sem, xlabel = "Short Poke Duration", xticks = 0:0.1:1, title = "Filtered: Poke dur. > 0.3, Interpoke dur. > 0.1")),
    (@df shortinter[2] bar(:Bin,:Count_mean, yerr = :Count_sem, xlabel = "Short Interpoke Duration")),
    layout = (2,1),
    size=(620,874)
    )
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusSelectionDurations.pdf"))
######### Filtering data
fCas_p = filter(r->r.PokeDur > 0.3 &&
    (r.Reward || ismissing(r.PostInterpoke) || (r.PostInterpoke > 0.1)) &&
    (r.Reward || r.PreInterpoke == 0 || ismissing(r.PreInterpoke) || (r.PreInterpoke > 0.1)),
    Cas_p)
    gd = groupby(fCas_p,[:MouseID,:Session,:Virus])
    fCas_s = combine(gd) do dd
        process_streaks(dd)
    end
    fpdf = group_frequency(fCas_s,:AfterLast, group = :Virus)
    # open_html_table(Cas_p)
    sort!(fpdf,:Xaxis)
    fpdf[isnan.(fpdf.Sem),:Sem] .= 0
    @df fpdf plot(:Xaxis, :Mean, ribbon = :Sem, group = :Virus,
        linecolor = :auto, ylims=(0,0.4), xlims = (0,29))
    pdf = group_frequency(Cas_s,:AfterLast, group = :Virus)
        sort!(pdf,:Xaxis)
        maximum(pdf.Xaxis)
        pdf[isnan.(pdf.Sem),:Sem] .= 0
        @df pdf plot(:Xaxis, :Mean, ribbon = :Sem, group = :Virus,
            linecolor = :auto, ylims=(0,0.4), xlims = (0,29))
plot(
    (@df pdf plot(:Xaxis, :Mean, ribbon = :Sem, group = :Virus,
        linecolor = :auto, ylims=(0,0.35), xlims = (0,29), xlabel = "Consecutive Failures",
        title = "Unfiltered")),
    (@df fpdf plot(:Xaxis, :Mean, ribbon = :Sem, group = :Virus,
        linecolor = :auto, ylims=(0,0.37), xlims = (0,29), xlabel = "Consecutive Failures",
        title = "Filtered: Poke dur. > 0.3, Interpoke dur. > 0.1")),
    layout = (2,1),
    size=(620,874)
    )
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusFilteredFailures.pdf"))
##
PALCas0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),fCas_s,Poisson())
PALCas1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),fCas_s,Poisson())
fCas_s.Check = predict(PALCas1)
open_html_table(fCas_s[:,[:MouseID,:Streak,:Virus,:AfterLast,:Check]])
Likelyhood_Ratio_test(PALCas0,PALCas1)
##
BCorrCas0 = fit(MixedModel,@formula(CorrectLeave ~ 1 + Streak + (1|MouseID)),fCas_s,Bernoulli())
BCorrCas1 = fit(MixedModel,@formula(CorrectLeave ~ 1 + Streak + Virus + (1|MouseID)),fCas_s,Bernoulli())
Likelyhood_Ratio_test(BCorrCas0,BCorrCas1)
##
#= example biting CD17
    trial 55 good clean minute x:xx
    trial 99 bad clean minute xx:xx
=#
streak_num = 55
if streak_num == 55
    pathstring = joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","CD17-GoodBiting.pdf")
    plttitle = "Good filtering"
elseif streak_num == 99
    pathstring = joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","CD17-BadBiting.pdf")
    plttitle = "Bad filtering"
end
tt = filter(r -> r.MouseID == "CD17" &&  r.Streak == streak_num, Cas_p)
    pltunfilt = plot(legend = false, xlabel = "Trials",
        yaxis = false, yticks = false, ylabel = "Time (seconds)",size=(800,800),
        title = "Unfiltered trial")
        for r in eachrow(tt)
            FLPDevelopment.session_plot!(pltunfilt, r)
        end
tt = filter(r -> r.MouseID == "CD17" &&  r.Streak == streak_num, fCas_p)
    pltfilt = plot(legend = false, xlabel = "Trials",
        yaxis = false, yticks = false, ylabel = "Time (seconds)",size=(800,800),
        title = "Filtered trial")
        for r in eachrow(tt)
            FLPDevelopment.session_plot!(pltfilt, r)
        end
plt = plot(
    pltunfilt,
    pltfilt,
    layout = (2,1),
    size=(620,874)
    )
savefig(plt,pathstring)
############## Poking time over trials
gd = groupby(Cas_p,[:Virus,:MouseID, :Streak])
df1 = combine(gd, :PokeDur => sum => :PokingTime, [:PokeIn,:PokeOut] => ((in,out) -> (out[end] - in[1])) => :TrialDuration)
df1[!,:PokingInTrial] = df1.PokingTime./df1.TrialDuration
filter!(r -> r.PokingInTrial< 1,df1)
gd1 = groupby(df1, [:Virus,:MouseID])
df2 = combine(gd1, :PokingInTrial => mean)
sort!(df2,:Virus)
open_html_table(df2)
gd2 = groupby(df2,:Virus)
df3 = combine(gd2, :PokingInTrial_mean => mean => :PokingInTrial_groupmean,
    :PokingInTrial_mean => sem => :PokingInTrial_groupsem)
open_html_table(df3)
############################### PSTH ########################

pokes_age = reallignpokes(Age_p)
@df filter(r -> r.PO_LR >= 0, pokes_age) ea_histogram(:PO_LR, bins = 20, xlims=(0,30))
psth_age = FLPDevelopment.pokes_psth(pokes_age.PI_LR,pokes_age.PO_LR)
bar(minimum(pokes_age.PI_LR):maximum(pokes_age.PO_LR) +1 ,psth_age,
    linewidth =0, color = :black,
    # xlims = (0,60),
    xticks = 0:300:4500,
    title = "Pokes PSTH Age-group",
    label = "All animals",
    xlabel = "Time from last reward poke out",
    ylabel = "Events count")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","PSTHAgepokesFull.pdf"))

pokes_cas = reallignpokes(Cas_p)
@df filter(r -> r.PO_LR >= 0, pokes_cas) ea_histogram(:PO_LR, bins = 20, xlims=(0,30))
psth_cas = FLPDevelopment.pokes_psth(pokes_cas.PI_LR,pokes_cas.PO_LR)
bar(minimum(pokes_cas.PI_LR):maximum(pokes_cas.PO_LR) +1 ,psth_cas,
    linewidth =0, color = :black,
    # xlims = (0,60),
    xticks = 0:300:4500,
    title = "Pokes PSTH Virus-group",
    label = "All animals",
    xlabel = "Time from last reward poke out",
    ylabel = "Events count")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","PSTHCaspokesFull.pdf"))
####### psth per animal
pokes_age = reallignpokes(Age_p)
gd = groupby(pokes_age, [:MouseID,:Age])
psth_cas = combine(x -> FLPDevelopment.pokes_psth2(x.PI_LR,x.PO_LR), gd)
for m in union(psth_cas.MouseID)
    dd = filter(r -> r.MouseID == m, psth_cas)
    g = dd[1,:Age]
    plt = @df dd bar(:Time.+1, :Psth,
        # xticks = 0:300:4500,
        xaxis=(:log10,[1,:auto]),
        linewidth =0,
        color = g == "Juveniles" ? get_color_palette(:auto,plot_color(:white))[2] : get_color_palette(:auto,plot_color(:white))[1],
        title = "Pokes PSTH Age-group",
        label = m,
        xlabel = "Time from last reward poke out",
        ylabel = "Events count")
        plt
    savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","PSTH",m*".png"))
end


pokes_cas = reallignpokes(Cas_p)
gd = groupby(pokes_cas, [:MouseID,:Virus])
psth_cas = combine(x -> FLPDevelopment.pokes_psth2(x.PI_LR,x.PO_LR), gd)
for m in union(psth_cas.MouseID)
    dd = filter(r -> r.MouseID == m, psth_cas)
    g = dd[1,:Virus]
    plt = @df dd bar(:Time.+1, :Psth,
        # xticks = 0:300:4500,
        xaxis=(:log10,[1,:auto]),
        linewidth =0,
        color = g == "Caspase" ? get_color_palette(:auto,plot_color(:white))[2] : get_color_palette(:auto,plot_color(:white))[1],
        title = "Pokes PSTH Virus-group",
        label = m,
        xlabel = "Time from last reward poke out",
        ylabel = "Events count")
        plt
    savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","PSTH",m*".png"))
end


####### psth per group
pokes_age = reallignpokes(Age_p)
gd = groupby(pokes_age, :Age)
psth_age = combine(gd,
    (:PI_LR,:PO_LR) => ((i,o) -> (FLPDevelopment.pokes_psth(i,o))) => :Psth,
    (:PI_LR,:PO_LR) => ((i,o) -> (collect(Int64(minimum(i)):1:Int64(maximum(o))))) => :time)
@df psth_age bar(:time, :Psth,
    group = :Age,
    xlims = (0,60),
    xticks = 0:10:60,
    linewidth =0, color = :auto,
    fillalpha = 0.5,
    title = "Pokes PSTH Age-group",
    xlabel = "Time from last reward poke out",
    ylabel = "Events count")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","GroupPSTHAge.pdf"))

pokes_cas = reallignpokes(Cas_p)
@df pokes_cas histogram(:PO_LR.+1, xaxis=(:log10,[1,:auto]))
gd = groupby(pokes_cas, :Virus)
psth_cas = combine(gd,
    (:PI_LR,:PO_LR) => ((i,o) -> (FLPDevelopment.pokes_psth(i,o))) => :Psth,
    (:PI_LR,:PO_LR) => ((i,o) -> (collect(Int64(minimum(i)):1:Int64(maximum(o))))) => :time)
@df psth_cas bar(:time.+ 1, :Psth,
    group = :Virus,
    # xlims = (0,60),
    xaxis=(:log10,[1,:auto]),
    xscale = :identity,
    xticks = :auto,
    linewidth = 0,
    color = :auto,
    fillalpha = 0.5,
    title = "Pokes PSTH Virus-group",
    xlabel = "Time from last reward poke out",
    ylabel = "Events count")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","GroupPSTHVirus.pdf"))

######## psth different bin size
pokes_cas = reallignpokes(Cas_p)
gd = groupby(pokes_cas, :MouseID)
psth_cas = combine(x -> FLPDevelopment.pokes_psth2(x.PI_LR,x.PO_LR), gd)
psth_cas = combine(x -> FLPDevelopment.pokes_psth2(x.PI_LR,x.PO_LR; bin_size = 0.5), gd)
psth_cas = combine(x -> FLPDevelopment.pokes_psth2(x.PI_LR,x.PO_LR; bin_size = 2 ), gd)



#=
fit psth
To fit the distribution a vector containing the time at which a poke was occurring.
For pokes that spans multiple time bins (seconds in this case) the value of each
bin has to be in the vector. To achieve that you need to collect from poke-in to
poke-out.
=#
test = vcat([collect(r.PI_LR:r.PO_LR) for r in eachrow(pokes_age)]...)
test = test[test .>= 0 .+ test.<=60] .+ 1
Gfit_age = fit(Gamma,test)
WMGfit_age = mixture_gamma_weighted(test)
WMEfit_age = mixture_exp_weighted(test)
WMEGfit_age = FLPDevelopment.mixture_gamma_exp_weighted(test)
Paretofit_age = fit(Pareto,test)


loglikelihood(Paretofit_age,test)
loglikelihood(WMGfit_age,test)
loglikelihood(WMEfit_age,test)
loglikelihood(WMEGfit_age,test)
loglikelihood(Gfit_age,test)

check = pdf.(Gfit_age,1:61)
check2 = pdf.(WMGfit_age,1:61)
check3 = pdf.(WMEfit_age,1:61)
check4 = pdf.(WMEGfit_age, 1:61)
check5 = pdf.(Paretofit_age, 1:61)
histogram(test)
plot!(0:60, check, label = "Gamma_fit", linecolor = :red, linewidth = 3)
plot!(0:60, check2, label = "Mixed Gamma_fit", linecolor = :cyan, linewidth = 3)
plot!(0:60, check3, label = "Mixed Exp_fit", linecolor = :green, linewidth = 3)
plot!(0:60, check3, label = "Mixed Gamma_Exp_fit", linecolor = :magenta, linewidth = 3)
plot!(0:60, check3, label = "Pareto_fit", linecolor = :yellow, linewidth = 3)
fit(Histogram,test)


savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","GammaPSTHAge.pdf"))
