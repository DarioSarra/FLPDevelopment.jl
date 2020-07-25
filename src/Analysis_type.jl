struct dvAnalysis
    mice_summary
    group_summary
    normality
    test
    plot
end

function dvAnalysis(df,xvar,yvar)
    if eltype(df[:,yvar]) == Bool
        println("Boolean vector, using fraction as summary")
        summary = fraction_true
    else
        summary = mean
    end
    df1 = mouse_summary(df,xvar,yvar; summary = summary)
    normality = test_normality(df1,xvar,yvar)
    df2 = group_summary(df1,xvar,yvar; normality = normality)
    test = test_difference(df1,xvar,yvar)
    plt = dvplot(df1,df2,xvar,yvar,test)
    dvAnalysis(df1, df2, normality, test, plt)
end
