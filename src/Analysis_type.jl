struct dvAnalysis
    mice_summary
    group_summary
    normality
    test
    plot
end

function dvAnalysis(df,xvar,yvar; yspan = :auto, ystep = :auto, summary_opt = :MEAN, bs = false, nonparametric = false)
    if eltype(df[:,yvar]) == Bool
        println("Boolean vector, using fraction as summary")
        summary_opt = :FRACTION_TRUE
    end
    if summary_opt == :MEAN
        summary = mean
    elseif summary_opt == :BOOTSTRAP_MEAN
        println("using bootstrap_mean as mouse summary")
        summary = bootstrap_mean
    elseif summary_opt == :FRACTION_TRUE
        summary = fraction_true
    elseif summary_opt == :SUM
        summary = +
    end
    df1 = mouse_summary(df,xvar,yvar; summary = summary)
    normality = test_normality(df1,xvar,yvar)
    if nonparametric
        norm = false
    else

    end
    df2 = group_summary(df1,xvar,yvar; normality = normality)
    test = test_difference(df1,xvar,yvar; normality = normality)
    plt = dvplot(df1,df2,xvar,yvar,test; yspan = yspan, ystep = ystep)
    dvAnalysis(df1, df2, normality, test, plt)
end

struct DoubleAnalysis
    mice_summary
    JarqueBera
    parametric_summary
    UnequalVarianceT
    parametric_plot
    nonparametric_summary
    MannWhitneyU
    nonparametric_plot
end

function DoubleAnalysis(df,xvar,yvar; yspan = :auto, ystep = :auto, summary_opt = :MEAN, bs = false)
    println("NEW APPROACH")
    if eltype(df[:,yvar]) == Bool
        println("Boolean vector, using fraction as summary")
        summary_opt = :FRACTION_TRUE
    end
    if summary_opt == :MEAN
        summary = mean
    elseif summary_opt == :BOOTSTRAP_MEAN
        println("using bootstrap_mean as mouse summary")
        summary = bootstrap_mean
    elseif summary_opt == :FRACTION_TRUE
        summary = fraction_true
    elseif summary_opt == :SUM
        summary = sum
    end
    df1 = mouse_summary(df,xvar,yvar; summary = summary)
    JarqueBera = test_normality(df1,xvar,yvar)
    parametric_summary = group_summary(df1,xvar,yvar; normality = true)
    UnequalVarianceT = test_difference(df1,xvar,yvar; normality = true)
    parametric_plot = dvplot(df1,parametric_summary,xvar,yvar,UnequalVarianceT; yspan = yspan, ystep = ystep)
    nonparametric_summary = group_summary(df1,xvar,yvar; normality = false)
    MannWhitneyU = test_difference(df1,xvar,yvar; normality = false)
    nonparametric_plot = dvplot(df1,nonparametric_summary,xvar,yvar,MannWhitneyU; yspan = yspan, ystep = ystep)

    DoubleAnalysis(df1, JarqueBera,
        parametric_summary, UnequalVarianceT, parametric_plot,
        nonparametric_summary, MannWhitneyU, nonparametric_plot)
end
