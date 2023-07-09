include(joinpath(splitpath(@__DIR__)[1:end-1]...,"Filters.jl"))
##############################
#= the possibility that individuals within the same group/condition are
using different strategies or are displaying different behavioral patterns
should be addressed=#
c = [x == "Adults" ? 1 : 2 for x in Age_s.Age]
    @df Age_s density(:LogDuration, group = :MouseID,
    linecolor = c,linealpha = 0.3,
    xlabel = "Leaving time (s)", ylabel = "Kernell density estimate",
    xticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), legend = false)
    c = [x == "Adults" ? 13 : 14 for x in Age_s.Age]
    @df Age_s density!(:LogDuration, group = :Age, linestyle = :dot,
    linewidth = 2, color = c, left_margin = -15px, bottom_margin = -15px)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_Bimodality.png"))
##
c = [x == "tdTomato" ? 1 : 2 for x in Cas_s.Virus]
    @df Cas_s density(:LogDuration, group = :MouseID,
    linecolor = c,linealpha = 0.3,
    xlabel = "Leaving time (s)", ylabel = "Kernell density estimate",
    xticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), legend = false)
    c = [x == "tdTomato" ? 13 : 14 for x in Cas_s.Virus]
    @df Cas_s density!(:LogDuration, group = :Virus, linestyle = :dot,
    linewidth = 2, color = c, left_margin = -15px, bottom_margin = -15px)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_Bimodality.png"))

##
Age_outliers = FLPDevelopment.outliersdf(Age_s, :LogDuration; group = :Age)
findall(Age_outliers.Outliers)
combine(groupby(Age_outliers,:Age)) do dd
    dc = countmap(length.(dd.Pos))
    DataFrame(N_peaks = collect(keys(dc)),N_mice = collect(values(dc)), Percent = collect(values(dc))./nrow(dd))
end
(12+11)/nrow(Age_outliers)
ecum = ecdf(Age_s.LogDuration)
r = range(0,3,length = 100)
x = 0.04
lowbin = Int64(ceil(x/Float64(r.step)))
highbin = lowbin +1
ecum(highbin*Float64(r.step)) - ecum(lowbin*Float64(r.step))


x(3)
Cas_outliers = FLPDevelopment.outliersdf(Cas_s, :LogDuration; group = :Virus)
findall(Cas_outliers.Outliers)
combine(groupby(Cas_outliers,:Virus)) do dd
    dc = countmap(length.(dd.Pos))
    DataFrame(N_peaks = collect(keys(dc)),N_mice = collect(values(dc)), Percent = collect(values(dc))./nrow(dd))
end
##
get_color_palette(:auto,17)[13:14]
@df Age_s violin(:MouseID, :LogDuration, side = :left, group = :Age,
        yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Leaving time density (s)",
        xticks = false, xlabel = "Mouse ID",
        size = (1500,900), ylims = (-1,3.5),
        color_palette = get_color_palette(:auto,17)[13:14]
)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_Violin.png"))
##
@df Cas_s violin(:MouseID, :LogDuration, side = :left, group = :Virus,
        yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Leaving time density (s)",
        xticks = false, xlabel = "Mouse ID",
        size = (1500,900), ylims = (-1,3.5),
        color_palette = get_color_palette(:auto,17)[13:14] 
)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_Violin.png"))

##
pdf_heatmap(Age_s,:Age,:LogDuration)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_Heatmap.png"))

pdf_heatmap(Cas_s,:Virus,:LogDuration)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_Heatmap.png"))
##
ddd = kde(Age_s.LogDuration)
sum(pdf(ddd,range(-1,3, length = 100)))
Age_ecdf = ecdf(Age_s.LogDuration)
Age_ecdf(2)
## KSampleAD test
cs = filter(r->r.Virus == "tdTomato",Cas_s)[:,[:MouseID,:LogDuration]]
coll = AbstractVector{<:Real}[]
combine(groupby(cs,:MouseID)) do dd
    push!(coll,dd.LogDuration)
end
coll
KSampleADTest(coll...)
a = [1,2,3]
v= [2,3,4]
KSampleADTest(xs::AbstractVector{<:Real}...; modified = true, nsim = 0)

#Show that mouse means distributions is not bimodal
Age_dd = combine(groupby(Age_s,[:Age,:MouseID]), :LogDuration => median => :LogDuration)
Age_AD = combine(groupby(Age_dd,[:Age])) do subdf
    dist = Normal(mean(subdf.LogDuration), std(subdf.LogDuration))
    ad = OneSampleADTest(subdf.LogDuration, dist)
    println(ad)
    (Obs = ad.n, Mean = ad.μ, SD = ad.σ, A= ad.A², Normal = pvalue(ad),)
end
Age_dd[!,:color] = [x=="Juveniles" ? 14 : 13 for x in Age_dd.Age]
sort!(Age_dd, :Age)
@df Age_dd dotplot(:Age, :LogDuration, markercolor = :color, legend = false, markersize = 5,
    yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Median Leaving time (s)", ylims = (0,2))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_LeavingTimeDist.png"))
##
Cas_dd = combine(groupby(Cas_s,[:Virus,:MouseID]), :LogDuration => mean => :LogDuration)
Cas_AD = combine(groupby(Cas_dd,[:Virus])) do subdf
    dist = Normal(mean(subdf.LogDuration), std(subdf.LogDuration))
    ad = OneSampleADTest(subdf.LogDuration, dist)
    println(ad)
    (Obs = ad.n, Mean = ad.μ, SD = ad.σ, A= ad.A², Normal = pvalue(ad),)
end
Cas_dd[!,:color] = [x=="Caspase" ? 14 : 13 for x in Cas_dd.Virus]
sort!(Cas_dd, :Virus)
@df Cas_dd dotplot(:Virus, :LogDuration, markercolor = :color, legend = false,markersize = 5,
    yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Median Leaving time (s)", ylims = (0,2))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_LeavingTimeDist.png"))
##
@df Age_dd boxplot(:Age, :LogDuration, legend = false, color = :lightgrey)
@df Age_dd dotplot!(:Age, :LogDuration, markercolor = :color, legend = false, markersize = 5,
yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Median Leaving time (s)", ylims = (0,2))
##
@df Cas_dd boxplot(:Virus, :LogDuration, legend = false, color = :lightgrey)
@df Cas_dd dotplot!(:Virus, :LogDuration, markercolor = :color, legend = false,markersize = 5,
    yticks = ([-1,0,1,2,3], ["0", "1", "10", "100", "1000"]), ylabel = "Median Leaving time (s)", ylims = (0,2))
##

