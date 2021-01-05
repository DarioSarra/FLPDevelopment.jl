function mouse_summary(df,xvar,yvar; summary = mean)
    gdc = groupby(df,[xvar,:MouseID])
    df1 = combine(yvar => summary => yvar,gdc)
    sort!(df1,xvar)
    firstval = union(df1[:,xvar])[1]
    df1[!,:xpos] = [v == firstval ? 1 : 2  for v in df1[:,xvar]]
    return df1
end
"""
    test_normality(df1,xvar,yvar)

Split a the vector yvar according to a binomial xvar vector in a DataFrame df1
Compute the Jarque-Bera statistic to test the null hypothesis that
both resulting vectors are normally distributed.
Return true if both are normally distributed.
"""
function test_normality(df1,xvar,yvar)
    cases = union(df1[:,xvar])
    case1 = df1[df1[:,xvar] .== cases[1], yvar]
    case2 = df1[df1[:,xvar] .== cases[2], yvar]
    n1 = pvalue(JarqueBeraTest(case1)) >= 0.05
    n2 = pvalue(JarqueBeraTest(case2)) >= 0.05
    # n1 = pvalue(ExactOneSampleKSTest(case1,Normal(mean(case1),std(case1)))) >= 0.05
    # n2 = pvalue(ExactOneSampleKSTest(case2,Normal(mean(case2),std(case2)))) >= 0.05
    all((n1,n2))
end

function group_summary(df1,xvar,yvar; normality = true)
    if normality
        central = mean
        confidence = OneSampleTTest
    else
        central = median
        confidence = SignedRankTest
    end
    df2 = combine(groupby(df1,xvar)) do dd
        m = round(central(dd[:,yvar]), digits = 2)
        ci = confint(confidence(dd[:,yvar]))
        ci1 = round(m - ci[1], digits = 2)
        ci2 = round(ci[2] - m, digits = 2)
        (Central = m, ERR = (ci1,ci2))
    end
    return df2
end

function test_difference(df1,xvar,yvar;normality = true)
    cases = union(df1[:,xvar])
    case1 = df1[df1[:,xvar] .== cases[1], yvar]
    case2 = df1[df1[:,xvar] .== cases[2], yvar]
    if normality
        test = UnequalVarianceTTest(case1,case2)
    else
        test = MannWhitneyUTest(case1,case2)
    end
    return test
end
"""
    `dvplot(df1,df2,xvar,yvar,test; yspan = :auto, ystep = :auto)`
df1 is the mouse summary dataframe
df2 is the group summary dataframe
xvar is the Symbol of the column to use on the x axes
yvar is the Symbol of the column to use on the y axes
"""
function dvplot(df1,df2,xvar,yvar,test; yspan = :auto, ystep = :auto)
    plt = @df df2 scatter(1:nrow(df2),:Central, yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),cols(xvar)),
        legend = false)
    @df df1 scatter!(:xpos,cols(yvar),
        markersize = 3,
        alpha = 0.5,
        color = :grey,
        ylims = yspan,
        yticks = ystep)
    if pvalue(test) < 0.05
        p = pvalue(test)
        message = p < 0.01 ? "p < 0.01" : "p < 0.05"
        if yspan == :auto
            m1 = maximum(df1[:,yvar])
            m2 = minimum(df1[:,yvar])
            span = (m1-m2)/20
            ref = m1 + span
        else
            span = (yspan[2] - yspan[1])/20
            ref = yspan[2] - span
        end

        plot!([1,1],[ref,ref+span])
        plot!([2,2],[ref,ref+span])
        plot!([1,2],[ref+span/2,ref+span/2])
        annotate!([(1.5,1.2*span+ref,
            Plots.text(message,
            7, :center))])
    end
    return plt
end

function check_cd9(df,xvar)
    p = plot(; label = false)
    for m in union(df.MouseID)
        lc = m == "CD09" ? :red : :grey
        @df filter(r -> r.MouseID == m,df) density!(cols(xvar), linecolor = lc,label = false)
    end
    p
end

function check_cd9(df,xvar, yvar)
    p = plot(; label = false)
    for m in union(df.MouseID)
        lc = m == "CD09" ? :red : :grey
        @df filter(r -> r.MouseID == m,df) scatter!(cols(xvar), cols(yvar), markercolor = lc,
        label = false, markersize = 2)
    end
    p
end
