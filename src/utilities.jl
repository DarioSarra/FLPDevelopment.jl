function fraction_true(v::AbstractArray{Bool})
    sum(v)/length(v)
end

function bootstrap_mean(v; sampling = BasicSampling, nsample = 1000)
    bs = bootstrap(mean, v, sampling(nsample))
    first(bs.t0)
end

function mouse_summary(df,xvar,yvar; summary = mean)
    gdc = groupby(df,[xvar,:MouseID])
    df1 = combine(yvar => summary => yvar,gdc)
    sort!(df1,xvar)
    firstval = union(df1[:,xvar])[1]
    df1[!,:xpos] = [v == firstval ? 1 : 2  for v in df1[:,xvar]]
    return df1
end

function test_normality(df1,xvar,yvar)
    cases = union(df1[:,xvar])
    case1 = df1[df1[:,xvar] .== cases[1], yvar]
    case2 = df1[df1[:,xvar] .== cases[2], yvar]
    n1 = pvalue(ExactOneSampleKSTest(case1,Normal(mean(case1),std(case1)))) >= 0.05
    n2 = pvalue(ExactOneSampleKSTest(case2,Normal(mean(case2),std(case2)))) >= 0.05
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
        m = central(dd[:,yvar])
        ci = confint(confidence(dd[:,yvar]))
        ci1 = m - ci[1]
        ci2 = ci[2] - m
        (Central = central(dd[:,yvar]), ERR = (ci1,ci2))
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

function dvplot(df1,df2,xvar,yvar,test; yspan = :auto)
    plt = @df df2 scatter(1:nrow(df2),:Central, yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),cols(xvar)),
        legend = false)
    @df df1 scatter!(:xpos,cols(yvar),
        markersize = 3,
        alpha = 0.5,
        color = :grey,
        ylims = yspan)
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
