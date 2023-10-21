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
Age_kde_mice = combine(groupby(Age_s, [:MouseID, :Age])) do dd
    k = kde(dd.LogDuration)
    (LeavingTime = collect(k.x), KernelDensity = k.density)
end
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","SFig3","A_Age_IndividualDensity.csv"), Age_kde_mice)
Age_kde = combine(groupby(Age_s, [:Age])) do dd
    k = kde(dd.LogDuration)
    (LeavingTime = collect(k.x), KernelDensity = k.density)
end
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","SFig3","A_Age_GrouplDensity.csv"), Age_kde)
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

Cas_kde_mice = combine(groupby(Cas_s, [:MouseID, :Virus])) do dd
    k = kde(dd.LogDuration)
    (LeavingTime = collect(k.x), KernelDensity = k.density)
end
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","SFig3","A_Cas_IndividualDensity.csv"), Cas_kde_mice)
Cas_kde = combine(groupby(Cas_s, [:Virus])) do dd
    k = kde(dd.LogDuration)
    (LeavingTime = collect(k.x), KernelDensity = k.density)
end
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","SFig3","A_Cas_GrouplDensity.csv"), Cas_kde)
## Animals PLeaving heatmaps
pdf_heatmap(Age_s,:Age,:LogDuration)
# savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_Heatmap.png"))

pdf_heatmap(Cas_s,:Virus,:LogDuration)
# savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_Heatmap.png"))