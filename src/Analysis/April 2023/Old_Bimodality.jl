include(joinpath(splitpath(@__DIR__)[1:end-1]...,"Filters.jl"))
##############################
#= the possibility that individuals within the same group/condition are
using different strategies or are displaying different behavioral patterns
should be addressed=#
## Animals KDE
c = [x == "Adults" ? 1 : 2 for x in Age_s.Age]
    @df Age_s density(:LogDuration, group = :MouseID,
    linecolor = c,linealpha = 0.3,
    xlabel = "Leaving time (s)", ylabel = "Kernell density estimate",
    xticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), legend = false, title = "Age manipulation")
    c = [x == "Adults" ? 13 : 14 for x in Age_s.Age]
    @df Age_s density!(:LogDuration, group = :Age, linestyle = :dot,
    linewidth = 2, color = c, left_margin = -15px, bottom_margin = -15px,
    )
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_Bimodality.png"))
##
c = [x == "tdTomato" ? 1 : 2 for x in Cas_s.Virus]
    @df Cas_s density(:LogDuration, group = :MouseID,
    linecolor = c,linealpha = 0.3,
    xlabel = "Leaving time (s)", ylabel = "Kernell density estimate",
    xticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), legend = false, title = "Virus manipulation")
    c = [x == "tdTomato" ? 13 : 14 for x in Cas_s.Virus]
    @df Cas_s density!(:LogDuration, group = :Virus, linestyle = :dot,
    linewidth = 2, color = c, left_margin = -15px, bottom_margin = -15px,
    )
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_Bimodality.png"))

## Animals PLeaving heatmaps
pdf_heatmap(Age_s,:Age,:LogDuration)
# savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_Heatmap.png"))

pdf_heatmap(Cas_s,:Virus,:LogDuration)
# savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_Heatmap.png"))


# One sample Anderson-Darling test. Show that mouse means distributions is not bimodal
Age_dd = combine(groupby(Age_s,[:Age,:MouseID]), :LogDuration => median => :LogDuration)
Age_AD = combine(groupby(Age_dd,[:Age])) do subdf
    dist = Normal(mean(subdf.LogDuration), std(subdf.LogDuration))
    ad = OneSampleADTest(subdf.LogDuration, dist)
    p_res = pvalue(ad)
    println(ad)
    (Obs = ad.n, Mean = ad.μ, SD = ad.σ, A= ad.A², P = p_res, Normal = p_res < 0.8,)
end
# open_html_table(Age_AD)
Age_dd[!,:color] = [x=="Juveniles" ? 14 : 13 for x in Age_dd.Age];
sort!(Age_dd, :Age)
@df Age_dd dotplot(:Age, :LogDuration, markercolor = :color, legend = false, markersize = 5,
    yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Median Leaving time (s)", ylims = (0,2))
# savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_LeavingTimeDist.png"))
Age_dd
##
Cas_dd = combine(groupby(Cas_s,[:Virus,:MouseID]), :LogDuration => median => :LogDuration)
Cas_AD = combine(groupby(Cas_dd,[:Virus])) do subdf
    dist = Normal(mean(subdf.LogDuration), std(subdf.LogDuration))
    ad = OneSampleADTest(subdf.LogDuration, dist)
    p_res = pvalue(ad)
    println(ad)
    (Obs = ad.n, Mean = ad.μ, SD = ad.σ, A= ad.A², P = p_res, Normal = p_res < 0.8,)
end
# open_html_table(Cas_AD)
Cas_dd[!,:color] = [x=="Caspase" ? 14 : 13 for x in Cas_dd.Virus];
sort!(Cas_dd, :Virus)
@df Cas_dd dotplot(:Virus, :LogDuration, markercolor = :color, legend = false,markersize = 5,
    yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Median Leaving time (s)", ylims = (0,2))
# savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_LeavingTimeDist.png"))

## One sample Anderson-Darling test by mouse
Age_mouse_AD = combine(groupby(Age_s,[:Age,:MouseID])) do subdf
    dist = Normal(mean(subdf.LogDuration), std(subdf.LogDuration))
    ad = OneSampleADTest(subdf.LogDuration, dist)
    p_res = pvalue(ad)
    println(ad)
    (Obs = ad.n, Mean = ad.μ, SD = ad.σ, A= ad.A², P = p_res, Normal = p_res < 0.8,)
end
open_html_table(Age_mouse_AD)

Cas_mouse_AD = combine(groupby(Cas_s,[:Virus,:MouseID])) do subdf
    dist = Normal(mean(subdf.LogDuration), std(subdf.LogDuration))
    ad = OneSampleADTest(subdf.LogDuration, dist)
    p_res = pvalue(ad)
    println(ad)
    (Obs = ad.n, Mean = ad.μ, SD = ad.σ, A= ad.A², P = p_res, Normal = p_res < 0.8,)
end
open_html_table(Cas_mouse_AD)

##
Age_outliers = FLPDevelopment.outliersdf(Age_s, :LogDuration; group = :Age)
findall(Age_outliers.Outliers)
Age_peaks = combine(groupby(Age_outliers,:Age)) do dd
    dc = countmap(length.(dd.Pos))
    DataFrame(N_peaks = collect(keys(dc)),N_mice = collect(values(dc)), Percent = collect(values(dc))./nrow(dd))
end

Cas_outliers = FLPDevelopment.outliersdf(Cas_s, :LogDuration; group = :Virus)
findall(Cas_outliers.Outliers)
Cas_peaks = combine(groupby(Cas_outliers,:Virus)) do dd
    dc = countmap(length.(dd.Pos))
    DataFrame(N_peaks = collect(keys(dc)),N_mice = collect(values(dc)), Percent = collect(values(dc))./nrow(dd))
end
##