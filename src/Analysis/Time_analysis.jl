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
#### Time Log Transformed
#### Poking after last reward
pokes_age = reallignpokes(Age_p)
filter!(r -> r.PO_LR > 0 && r.PI_LR > 0, pokes_age)
transform!(pokes_age, :PI_LR => (x -> log10.(x)) => :Log_PI)
transform!(pokes_age, :PO_LR => (x -> log10.(x)) => :Log_PO)
gd = groupby(pokes_age, [:Age,:MouseID])
psth_age = combine(x -> FLPDevelopment.pokes_psth(x.Log_PI,x.Log_PO; bin_size = 0.1), gd)

for m in union(psth_age.MouseID)
    dd = filter(r -> r.MouseID == m, psth_age)
    g = dd[1,:Age]
    plt = @df dd bar(:Time, :Psth,
        # xticks = (1:9:length(Bas),Bas[1:9:length(Bas)]),
        # xaxis=(:log10,[0.1,:auto]),
        linewidth =1,
        color = g == "Adults" ? get_color_palette(:auto,plot_color(:white))[2] : get_color_palette(:auto,plot_color(:white))[1],
        title = "Age-group",
        label = m,
        xlabel = "Time (logscale s)",
        ylabel = "Events count")
    savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","LogTransformPSTH",m*".png"))
end


pokes_cas = reallignpokes(Cas_p)
filter!(r -> r.PO_LR > 0 && r.PI_LR > 0, pokes_cas)
transform!(pokes_cas, :PI_LR => (x -> log10.(x)) => :Log_PI)
transform!(pokes_cas, :PO_LR => (x -> log10.(x)) => :Log_PO)
gd = groupby(pokes_cas, [:Virus,:MouseID])
psth_cas = combine(x -> FLPDevelopment.pokes_psth(x.Log_PI,x.Log_PO; bin_size = 0.1), gd)

for m in union(psth_cas.MouseID)
    dd = filter(r -> r.MouseID == m, psth_cas)
    g = dd[1,:Virus]
    plt = @df dd bar(:Time, :Psth,
        # xticks = (1:9:length(Bas),Bas[1:9:length(Bas)]),
        # xaxis=(:log10,[0.1,:auto]),
        linewidth =1,
        color = g == "Caspase" ? get_color_palette(:auto,plot_color(:white))[2] : get_color_palette(:auto,plot_color(:white))[1],
        title = "Virus-group",
        label = m,
        xlabel = "Time (logscale s)",
        ylabel = "Events count")
    savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","LogTransformPSTH",m*".png"))
end

### Log Interpoke
pokes_age = reallignpokes(Age_p)
mplt_age = FLPDevelopment.bymouse_logtimehist(pokes_age,:PostInterpoke;digits = 1)
pokes_cas = reallignpokes(Cas_p)
mplt_cas = FLPDevelopment.bymouse_logtimehist(pokes_cas,:PostInterpoke;digits = 1)
filter!(r -> r.PO_LR > 0 && r.PI_LR > 0 &&
    !ismissing(r.PostInterpoke), pokes_age)
transform!(pokes_age, :PostInterpoke => (x -> round.(log10.(x), digits = 1)) => :Log_PostInterpoke)
pokes_age[findall(pokes_age.Log_PostInterpoke .== 0),:Log_PostInterpoke] .= 0.0
# open_html_table(sort(pokes_age,:Log_PostInterpoke))
gd = groupby(pokes_age, [:Age,:MouseID, :Log_PostInterpoke])
interpoke_age = combine(gd, nrow => :Count)
union(interpoke_age.Age)

for m in union(interpoke_age.MouseID)
    dd = filter(r -> r.MouseID == m, interpoke_age)
    g = dd[1,:Age]
    sort!(dd,:Log_PostInterpoke)
    plt = @df dd bar(:Log_PostInterpoke, :Count,
        linewidth =1,
        xlims = (-3,3),
        color = g == "Juveniles" ? get_color_palette(:auto,plot_color(:white))[2] : get_color_palette(:auto,plot_color(:white))[1],
        title = "Age-group",
        label = m,
        xlabel = "Interpoke duration (logscale s)",
        ylabel = "Events count")
    savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","Interpoke",m*".png"))
end

pokes_cas = reallignpokes(Cas_p)
filter!(r -> r.PO_LR > 0 && r.PI_LR > 0 &&
    !ismissing(r.PostInterpoke), pokes_cas)
transform!(pokes_cas, :PostInterpoke => (x -> round.(log10.(x), digits = 1)) => :Log_PostInterpoke)
pokes_cas[findall(pokes_cas.Log_PostInterpoke .== 0),:Log_PostInterpoke] .= 0.0
gd = groupby(pokes_cas, [:Virus,:MouseID, :Log_PostInterpoke])
interpoke_cas = combine(gd, nrow => :Count)
sort!(interpoke_cas,[:MouseID,:Log_PostInterpoke])
union(interpoke_cas.Virus)
for m in union(interpoke_cas.MouseID)
    dd = filter(r -> r.MouseID == m, interpoke_cas)
    g = dd[1,:Virus]
    plt = @df dd bar(:Log_PostInterpoke, :Count,
        linewidth =1,
        xlims = (-3,3),
        color = g == "Caspase" ? get_color_palette(:auto,plot_color(:white))[2] : get_color_palette(:auto,plot_color(:white))[1],
        title = "Virus-group",
        label = m,
        xlabel = "Interpoke duration (logscale s)",
        ylabel = "Events count")
    savefig(plt,joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","Interpoke",m*".png"))
end
############ Interpoke group
pokes_age = reallignpokes(Age_p)
filter!(r -> r.PO_LR > 0 && r.PI_LR > 0 &&
    !ismissing(r.PostInterpoke), pokes_age)
transform!(pokes_age, :PostInterpoke => (x -> round.(log10.(x), digits = 1)) => :Log_PostInterpoke)
pokes_age[findall(pokes_age.Log_PostInterpoke .== 0),:Log_PostInterpoke] .= 0.0
gd = groupby(pokes_age, [:Age,:MouseID, :Log_PostInterpoke])
interpoke_age = combine(gd, nrow => :Count)
interpokegroup_age = combine(groupby(interpoke_age,[:Age,:MouseID]), :Log_PostInterpoke, :Count => (x -> x./sum(x)) => :Fraction)
sort!(interpokegroup_age,[:Age,:MouseID, :Log_PostInterpoke])
interpokegroup_age = combine(groupby(interpokegroup_age,[:Age,:Log_PostInterpoke]), :Fraction => mean, :Fraction => sem)
sort!(interpokegroup_age,[:Age, :Log_PostInterpoke])
filter!(r -> !isnan(r.Fraction_sem), interpokegroup_age)
@df interpokegroup_age bar(:Log_PostInterpoke , :Fraction_mean, group = :Age, yerr = :Fraction_sem, alpha = 0.5)

pokes_cas = reallignpokes(Cas_p)
filter!(r -> r.PO_LR > 0 && r.PI_LR > 0 &&
    !ismissing(r.PostInterpoke), pokes_cas)
transform!(pokes_cas, :PostInterpoke => (x -> round.(log10.(x), digits = 1)) => :Log_PostInterpoke)
pokes_cas[findall(pokes_cas.Log_PostInterpoke .== 0),:Log_PostInterpoke] .= 0.0
gd = groupby(pokes_cas, [:Virus,:MouseID, :Log_PostInterpoke])
interpoke_cas = combine(gd, nrow => :Count)
sort!(interpoke_cas,[:MouseID,:Log_PostInterpoke])
interpokegroup_cas = combine(groupby(interpoke_cas,[:Virus,:MouseID]), :Log_PostInterpoke, :Count => (x -> x./sum(x)) => :Fraction)
sort!(interpokegroup_cas,[:Virus,:MouseID, :Log_PostInterpoke])
interpokegroup_cas = combine(groupby(interpokegroup_cas,[:Virus,:Log_PostInterpoke]), :Fraction => mean, :Fraction => sem)
sort!(interpokegroup_cas,[:Virus, :Log_PostInterpoke])
filter!(r -> !isnan(r.Fraction_sem), interpokegroup_cas)
@df interpokegroup_cas bar(:Log_PostInterpoke , :Fraction_mean, group = :Virus, yerr = :Fraction_sem, alpha = 0.5)
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
############ example biting CD19 minute 3:25
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
###
FLPDevelopment.all_duration_analysis(Age_p, Age_s, :Age)[1]
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","AgeTrial.pdf"))
FLPDevelopment.all_duration_analysis(Age_p, Age_s, :Age)[2]
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","AgePoke.pdf"))
FLPDevelopment.all_duration_analysis(Cas_p, Cas_s, :Virus)[1]
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","CasTrial.pdf"))
FLPDevelopment.all_duration_analysis(Cas_p, Cas_s, :Virus)[2]
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","DistTime","CasPoke.pdf"))
