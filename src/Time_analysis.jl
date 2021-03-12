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
plt = plot(trial..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Trial Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusTrialDurations.pdf"))

check = filter(r -> r.Travel_to > 0, Cas_s)
travel = FLPDevelopment.bin_duration(check,:Travel_to,:Virus)
plt = plot(travel..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Travel Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusTravelDurations.pdf"))

check = filter(r -> r.PreInterpoke > 0, Cas_p)
inter = FLPDevelopment.bin_duration(check,:PreInterpoke,:Virus)
check2 = FLPDevelopment.calculate_bin_duration(check,:PreInterpoke,:Virus)
bn = round.(union(check2.Bin), digits =2)
bn = join(string.(bn), " ")
plt = plot(inter..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "InterPoke Duration\n" * bn)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusInterpokeDurations.pdf"))

check = filter(r -> 0 < r.PreInterpoke <= 1, Cas_p)
inter = FLPDevelopment.bin_duration(check,:PreInterpoke,:Virus; modality = :STEP, unit_step = 0.05)
plt = plot(inter..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "InterPoke Duration\n" * bn)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusShortInterpokeDurations.pdf"))

check = filter(r -> r.Reward, Cas_p)
rew = FLPDevelopment.bin_duration(check,:PokeDur,:Virus)
plt = plot(rew..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Rewards Poke Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusRewardDurations.pdf"))

check = filter(r -> !r.Reward, Cas_p)
failures = FLPDevelopment.bin_duration(check,:PokeDur,:Virus)
plt = plot(failures..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Failures Poke Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusFailuresDurations.pdf"))

check = filter(r -> 0 < !r.Reward && r.PokeDur <= 0.5, Cas_p)
failures = FLPDevelopment.bin_duration(check,:PokeDur,:Virus; modality = :STEP, unit_step = 0.05)
plt = plot(failures..., size = (1754,1240), layout = (3,5))
maintitle!(plt, "Short Failures Poke Duration")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","VirusShortFailuresDurations.pdf"))
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
#########
fCas_p = filter(r->r.PokeDur > 0.35 &&
    (ismissing(r.PostInterpoke) || (r.PostInterpoke > 0.1)) &&
    (r.PreInterpoke == 0 || ismissing(r.PreInterpoke) || (r.PreInterpoke > 0.1)),
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
##
pdf = group_frequency(Cas_s,:AfterLast, group = :Virus)
sort!(pdf,:Xaxis)
maximum(pdf.Xaxis)
pdf[isnan.(pdf.Sem),:Sem] .= 0
@df pdf plot(:Xaxis, :Mean, ribbon = :Sem, group = :Virus,
    linecolor = :auto, ylims=(0,0.4), xlims = (0,29))
##
ALCas0 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak  + (1|MouseID)),fCas_s))
ALCas1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),fCas_s))
Likelyhood_Ratio_test(ALCas0,ALCas1)
PALCas0 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + (1|MouseID)),fCas_s,Poisson())
PALCas1 = fit(MixedModel,@formula(AfterLast ~ 1 + Streak + Virus + (1|MouseID)),fCas_s,Poisson())
Likelyhood_Ratio_test(PALCas0,PALCas1)
##
# example biting CD17 minute x:xx
tt = filter(r -> r.MouseID == "CD17" &&  r.Streak == 99, Cas_p)
plt = plot(legend = false, xlabel = "Trials",
    yaxis = false, yticks = false, ylabel = "Time (seconds)",size=(800,800))
    # xlims = (38,40), xticks = 0:0.2:300)
    for r in eachrow(tt)
        FLPDevelopment.session_plot!(plt, r)
    end
    plt
savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","CD17-biting.pdf"))
tt = filter(r -> r.MouseID == "CD17" &&  r.Streak == 99, fCas_p)
plt = plot(legend = false, xlabel = "Trials",
    yaxis = false, yticks = false, ylabel = "Time (seconds)",size=(800,800))
    # xlims = (38,40), xticks = 0:0.2:300)
    for r in eachrow(tt)
        FLPDevelopment.session_plot!(plt, r)
    end
    plt
savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","CD17-bitingcleaned.pdf"))
##
check = filter(r -> 0 < !r.Reward && r.PokeDur <= 1, Cas_p)
failures = FLPDevelopment.calculate_bin_duration(check,:PokeDur,:Virus; modality = :STEP, unit_step = 0.05)
sort!(failures,[:Bin,:MouseID])
groupdf = combine(groupby(failures,:Bin), :Count=> mean, :Count => sem)
@df groupdf bar(:Bin, :Count_mean, yerror = :Count_sem)
