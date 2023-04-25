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
Cas_outliers = FLPDevelopment.outliersdf(Cas_s, :LogDuration; group = :Virus)
findall(Cas_outliers.Outliers)
open_html_table(Cas_outliers)
##
function outliersdf(df,var;group = nothing)
    mouse = union(df.MouseID)[2]
        v = filter(r -> r.MouseID == mouse, df)[:,var]
        d = filter(r -> r.MouseID != mouse, df)[:,var]
        pvalue(ApproximateTwoSampleKSTest(v, d))
    isnothing(group) ? (df1 = df) : (df1 = groupby(df,[group]))
    test = combine(df1) do dd
        res = Float64[]
        mouse = String[]
        for m in unique(dd.MouseID)
            push!(mouse,m)
            d = filter(r -> r.MouseID != m, dd)[:,var]
            v = filter(r -> r.MouseID == m, df)[:,var]
            push!(res,pvalue(ApproximateTwoSampleKSTest(v, d)))
        end
        return (MouseID = mouse, P = res)
    end
    transform!(test, :P => ByRow(x-> x < 0.05) => :Outlier)
    return test
end

Age_test = outliersdf(Age_s,:LogDuration; group = :Age)
Age_outliers = filter(r-> r.Outlier == 1, Age_test).MouseID
Cas_test = outliersdf(Cas_s,:LogDuration; group = :Virus)
Cas_outliers = filter(r-> r.Outlier == 1, Cas_test).MouseID
##
