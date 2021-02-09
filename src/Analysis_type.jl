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
    df1 = individual_summary(df,xvar,yvar; summary = summary)
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
    bplot
    vplot
end

function DoubleAnalysis(df,xvar,yvar; yspan = :auto, ystep = :auto, summary_opt = :MEAN, showmice = false)
    df1 = step1_DoubleAnalysis(df,xvar,yvar; summary_opt = summary_opt)
    step2_DoubleAnalysis(df1,xvar,yvar; yspan = yspan, ystep = ystep, summary_opt = summary_opt, showmice = false)
end

function step1_DoubleAnalysis(df,xvar,yvar; summary_opt = :MEAN)
    if eltype(df[:,yvar]) == Bool
        if yvar == :IncorrectLeave
            println("IncorrectLeave vector, using incorrect_fraction as summary")
            summary_opt = :FRACTION_INCORRECT
        else
            println("Boolean vector, using fraction as summary")
            summary_opt = :FRACTION_TRUE
        end
    end
    if summary_opt == :MEAN
        summary = mean
    elseif summary_opt == :BOOTSTRAP_MEAN
        println("using bootstrap_mean as mouse summary")
        summary = bootstrap_mean
    elseif summary_opt == :FRACTION_TRUE
        summary = fraction_true
    elseif summary_opt == :FRACTION_INCORRECT
        summary = incorrect_fraction
    elseif summary_opt == :SUM
        summary = sum
    elseif summary_opt == :MEDIAN
        summary = median
    elseif summary_opt == :MODE
        summary = mode
    end
    df1 = individual_summary(df,xvar,yvar; summary = summary)
end
function step2_DoubleAnalysis(df1,xvar,yvar; yspan = :auto, ystep = :auto, summary_opt = :MEAN, showmice = false)
    JarqueBera = test_normality(df1,xvar,yvar)
    parametric_summary = group_summary(df1,xvar,yvar; normality = true)
    UnequalVarianceT = test_difference(df1,xvar,yvar; normality = true)
    parametric_plot = dvplot(df1,parametric_summary,xvar,yvar,UnequalVarianceT; yspan = yspan, ystep = ystep, showmice = showmice)
    nonparametric_summary = group_summary(df1,xvar,yvar; normality = false)
    MannWhitneyU = test_difference(df1,xvar,yvar; normality = false)
    nonparametric_plot = dvplot(df1,nonparametric_summary,xvar,yvar,MannWhitneyU; yspan = yspan, ystep = ystep, showmice = showmice)
    flipaxis = xvar == :Virus
    bplot = @df df1 boxplot(string.(cols(xvar)),cols(yvar), legend = false, xflip = flipaxis)
    if length(union(df1[:,xvar])) != 2
        vplot = @df df1 violin(string.(cols(xvar)),cols(yvar), legend = false, xflip = flipaxis)
    else
        vals = union(df1[:,xvar])
        vplot = @df  df1[df1[:,xvar] .== vals[1],:] violin(string.(cols(xvar)),cols(yvar),
             label = vals[1], side = :left, legend = :top)
         @df df1[df1[:,xvar] .== vals[2],:] violin!(string.(cols(xvar)),cols(yvar),
            label = vals[2], side = :right)
        @df df1[df1[:,xvar] .== vals[1],:] dotplot!(string.(cols(xvar)),cols(yvar),
            label = "", side = :left, legend = :top, markersize = 4, markercolor = :grey, markeralpha = 0.7)
        @df df1[df1[:,xvar] .== vals[2],:] dotplot!(string.(cols(xvar)),cols(yvar),
            label = "", side = :right, legend = :top, markersize = 4, markercolor = :grey,
            markeralpha = 0.7)
    end
    DoubleAnalysis(df1, JarqueBera,
        parametric_summary, UnequalVarianceT, parametric_plot,
        nonparametric_summary, MannWhitneyU, nonparametric_plot,
        bplot, vplot)
end
