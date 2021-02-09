"""
    `summary_xy(df,xvar,yvar; summary = mean, err = :MouseID, group = nothing)`
"""
function summary_xy(df,xvar,yvar; summary = mean, err = :MouseID, group = nothing)
    if group == nothing
        m_res = individual_summary(df,xvar,yvar; summary = summary, err = err)
        gd = groupby(m_res,xvar)
        res = combine(gd,yvar => mean => :Mean, yvar => sem => :Sem)
        return dropnan!(res)
    else
        gd1 = groupby(df,group)
        res = combine(gd1) do dd
            m_res = individual_summary(dd,xvar,yvar; summary = summary, err = err)
            gd2 = groupby(m_res,xvar)
            combine(gd2,yvar => mean => :Mean, yvar => sem => :Sem)
        end
        return dropnan!(res)
    end
end

function individual_kde(df,var; err = :MouseID, points = 100, bounds = extrema(df[:,var]))
    # axis = kde(df[:,var], boundary = bounds, npoints = points).x
    axis = range(bounds...,length = points)
    gd1 = groupby(df, err)
    df1 = combine(gd1) do dd
        ka = kde(dd[:,var], npoints = points#=, boundary = bounds=#)
        (Xaxis = collect(axis), Vals = [pdf(ka,x) for x in axis])
    end
end

function group_kde(df,var; err = :MouseID, group = nothing, points = 100, bounds = extrema(df[:,var]))
    if isnothing(group)
        df = individual_kde(df,var; err = err, points = points, bounds = bounds)
        gd = groupby(df,:Xaxis)
        return combine(gd, :Vals => mean => :Mean, :Vals => sem => :Sem)
    else
        gd1 = groupby(df, group)
        res = combine(gd1) do dd
            m_res = individual_kde(dd,var; err = err, points = points, bounds = bounds)
            gd2 = groupby(m_res,:Xaxis)
            combine(gd2,:Vals => mean => :Mean, :Vals => sem => :Sem)
        end
    end
end

function individual_cdf(df,var; err = :MouseID, group = nothing, points = 100, bounds = extrema(df[:,var]))
    axis = ecdf(df[:,var]).sorted_values
    gd1 = groupby(df, err)
    df1 = combine(gd1) do dd
        CUM = ecdf(dd[:,var])
        (Xaxis = collect(axis), Vals = CUM(axis))
    end
end

function group_cdf(df,var; err = :MouseID, group = nothing, points = 100, bounds = extrema(df[:,var]))
    if isnothing(group)
        df = individual_cdf(df,var; err = err, points = points, bounds = bounds)
        gd = groupby(df,:Xaxis)
        return combine(gd, :Vals => mean => :Mean, :Vals => sem => :Sem)
    else
        gd1 = groupby(df, group)
        res = combine(gd1) do dd
            m_res = individual_cdf(dd,var; err = err, points = points, bounds = bounds)
            gd2 = groupby(m_res,:Xaxis)
            combine(gd2,:Vals => mean => :Mean, :Vals => sem => :Sem)
        end
    end
end

function individual_frequency(df,xvar; err = :MouseID)
    axis = minimum(df[:,xvar]):maximum(df[:,xvar])
    gd1 = groupby(df, err)
    df1 = combine(gd1) do dd
        df2 = DataFrame(Xaxis = Int[], Vals = Float64[])
        for (k,i) in countmap(dd[:,xvar])
            push!(df2, [k,i/nrow(dd)])
        end
        df2
    end
    return sort(df1, :Xaxis)
end

function group_frequency(df,var; err = :MouseID, group = nothing)
    if isnothing(group)
        df = individual_frequency(df,var; err = err)
        gd = groupby(df,:Xaxis)
        return combine(gd, :Vals => mean => :Mean, :Vals => sem => :Sem)
    else
        gd1 = groupby(df, group)
        res = combine(gd1) do dd
            m_res = individual_frequency(dd,var; err = err)
            gd2 = groupby(m_res,:Xaxis)
            combine(gd2,:Vals => mean => :Mean, :Vals => sem => :Sem)
        end
    end
end

group_distribution(df, var; args...) = length(unique(df[:,var])) > 40 ? group_kde(df, var; args ...) : group_frequency(df, var; args ...)

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
    # n2 = pvalue(ExactOneSampleKSTest(case2,Normal(mean(case2),std(case2)))) >= 0.05
    # n1 = pvalue(ExactOneSampleKSTest(case1,Normal(mean(case1),std(case1)))) >= 0.05
    # n1 = pvalue(OneSampleADTest(case1, Normal())) >= 0.05
    # n2 = pvalue(OneSampleADTest(case2,Normal())) >= 0.05
    all([n1,n2])
end

function individual_summary(df,xvar,yvar; summary = mean, err = :MouseID)
    gdc = groupby(df,[xvar,err])
    df1 = combine(yvar => summary => yvar,gdc)
    sort!(df1,xvar)
    firstval = union(df1[:,xvar])[1]
    df1[!,:xpos] = [v == firstval ? 1 : 2  for v in df1[:,xvar]]
    return df1
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
        # m = round(central(dd[:,yvar]), digits = 2)
        # ci = confint(confidence(dd[:,yvar]))
        # ci1 = round(m - ci[1], digits = 2)
        # ci2 = round(ci[2] - m, digits = 2)
        b = bootstrap(central,dd[:,yvar], BasicSampling(100000))
        m, lower, upper = confint(b,BasicConfInt(0.95))[1]
        ci1 = m - lower#, digits = 2
        ci2 = upper - m#, digits = 2
        (Central = m, ERR = (ci2,ci1))
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
df1 is the individual summary dataframe
df2 is the group summary dataframe
xvar is the Symbol of the column to use on the x axes
yvar is the Symbol of the column to use on the y axes
"""
function dvplot(df1,df2,xvar,yvar,test; yspan = :auto, ystep = :auto, showmice = false)
    plt = @df df2 scatter(1:nrow(df2),:Central,
        yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),cols(xvar)),
        legend = false)
    if showmice
        @df df1 scatter!(:xpos,cols(yvar),
        markersize = 3,
        alpha = 0.5,
        color = :grey,
        ylims = yspan,
        yticks = ystep)
    end
    if pvalue(test) < 0.05
        p = pvalue(test)
        message = p < 0.01 ? "p < 0.01" : "p < 0.05"
        if yspan == :auto
            if showmice #use mice maximum and minimum values to set the plot span
                m1 = maximum(df1[:,yvar])
                m2 = minimum(df1[:,yvar])
            else #use yerror minimu and maximum values to set the plot span
                m1 = maximum(df2.Central .+ last.(df2.ERR))
                m2 = maximum(df2.Central .- first.(df2.ERR))
            end
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
## check distributions
function check_distribution(df, var; grouping = nothing, summary_opt = :MEAN)
    if isnothing(grouping)
        if in(:Age,propertynames(df))
            grouping = :Age
        elseif in(:Virus,propertynames(df))
            grouping = :Virus
        else
            error("Enable to find grouping variable")
         end
    end
    # Density
    grouped_df = group_distribution(df,var, group = grouping)
    filter!(r -> !isnan(r.Sem), grouped_df)
    density_plot = @df grouped_df plot(:Xaxis,:Mean, group = cols(grouping), xlabel = string(var), ribbon = :Sem, linecolor = :auto)
    q95 = quantile(df[:,var],0.95)
    vline!([q95], line = :dash)
    allignment = q95 > median(density_plot[1][1][:x]) ? :right : :left
    value = q95 > 1 ? string(Int64(round(q95, digits = 0))) : string(round(q95, digits = 2))
    annotate!(q95, maximum(density_plot[1][1][:y])/2,Plots.text(" 95th percentile: "* value * " \n ndata: " * string(nrow(df)) * " ",10,allignment))
    # Value over trials
    df[!,:BinnedStreak] = bin_axis(df.Streak; unit_step = 5)
    dfSum = summary_xy(df,:BinnedStreak,var; summary = mean, err = :MouseID, group = grouping)
    overtrialplot = @df dfSum plot(:BinnedStreak,:Mean, ribbon = :Sem, group = cols(grouping), linecolor = :auto,
    legend = :topleft, ylabel = String(var), xlabel = "Trial")
    # test and violin plot
    an_res = DoubleAnalysis(df,grouping,var, summary_opt = summary_opt)
    density_plot, overtrialplot, an_res.nonparametric_plot, an_res.vplot
end

function check_distributions(df_s,df_p; grouping = nothing, summary_opt = :MEAN)
    if isnothing(grouping)
        if in(:Age,propertynames(df_s))
            grouping = :Age
        elseif in(:Virus,propertynames(df_s))
            grouping = :Virus
         end
    end
    # AfterLast
    dAF, tAF, sAF, bAF = check_distribution(df_s,:AfterLast)
    # Error
    dErr, tErr, sErr, bErr = check_distribution(df_s,:IncorrectLeave)
    # Trial Rewards Rate
    df_s[!,:BinnedStreak] = bin_axis(df_s.Streak; unit_step = 5)
    dfS = summary_xy(df_s,:BinnedStreak,:Cum_Rewards; summary = mean, err = :MouseID, group = grouping)
    tTRR = @df dfS plot(:BinnedStreak,:Mean, ribbon = :Sem, group = cols(grouping), linecolor = :auto,
    legend = :topleft, ylabel = "Amount of Rewards", xlabel = "Trial")
    gdc = groupby(df_s,[grouping,:MouseID])
    df1 = combine(gdc, [:Num_Rewards, :Streak] => ((r,s) -> sum(r)/maximum(s)) => :RewxTrial)
    D = step2_DoubleAnalysis(df1,grouping,:RewxTrial)
    sTRR, bTRR = D.nonparametric_plot, D.vplot
    # Travel
    dTV, tTV, sTV, bTV = check_distribution(df_s,:Travel_to; summary_opt = summary_opt)
    #Interpoke
    dIP, tIP,sIP, bIP = check_distribution(df_p,:PreInterpoke; summary_opt = summary_opt)
    #Duration
    dDT, tDT, sDT, bDT = check_distribution(df_s,:Trial_duration; summary_opt = summary_opt)

    p1 = plot(tAF, sAF, bAF,
        tErr, sErr, bErr,
        tTRR, sTRR, bTRR,
        layout = grid(3,3),
        size=(1200,1200),
        thickness_scaling = 1)
    p2 = plot(tDT, sDT, bDT,
        tIP, sIP, bIP,
        tTV, sTV, bTV,
        layout = grid(3,3),
        size=(1200,1200),
        thickness_scaling = 1)
        (p1,p2)
end

function maintitle!(plt,what::String)
    y = ones(3)
    title = Plots.scatter(y, marker=0, markeralpha=0, grid = false,
    axis=false, leg=false,ticks = nothing, size = (200,1200),
    annotations=(2, y[2], Plots.text(what)))
    plot(title, plt, layout = grid(2,1, heights = [0.05,0.95]))
end
