"""
    `individual_summary(df,xvar,yvar; summary = mean, err = :MouseID)`
    computes a summary statistic (default = mean) for each individual
    (default = :MouseID) used to later calculate the error across individual
"""

function individual_summary(df,xvar,yvar; summary = mean, err = :MouseID)
    #evaluate if it's going over multiple grouping columns
    multiple = typeof(xvar) <: AbstractVector && sizeof(xvar) > 1
    if multiple
        gdc = groupby(df,vcat(xvar,err))
    else
        gdc = groupby(df,[xvar,err])
    end
    df1 = combine(gdc,yvar => summary => yvar)
    sort!(df1,xvar)
    if multiple
        label_df = unique(df1[:,xvar])
    else
        label_df = unique(df1[:,[xvar]])
    end
    label_df[!,:xpos] = 1:nrow(label_df)
    leftjoin(df1,label_df; on = xvar)
    try
        firstval = union(df1[:,xvar])[1]
        df1[!,:xpos] = [v == firstval ? 1 : 2  for v in df1[:,xvar]]
    catch
        println("Can't define x position")
    end
    return df1
end

"""
    `group_summary(df1,xvar,yvar; normality = true)`
    using bootstrap computes central measure plus 95% CI
    if normality is true it uses mean as central measure
    else it uses median
"""

function group_summary(df1,xvar,yvar; normality = true)
    if normality
        central = mean
    else
        central = median
    end
    df2 = combine(groupby(df1,xvar)) do dd
        # using bootstrap to calculate the central point and the 95% CI
        b = bootstrap(central, dd[:,yvar], BasicSampling(100000))
        #bootstrap confint returns a tuple containing a tuple with the statistic
        #estimates and the CI lower and uppper bounds
        m, lower, upper = confint(b,PercentileConfInt(0.95))[1]
        ci1 = m - lower
        ci2 = upper - m
        (Central = m, ERR = (ci1,ci2))
    end
    return df2
end

function MWU_test(df1,group,yvar)
    cases = union(df1[:,group])
    case1 = df1[df1[:,group] .== cases[1], yvar]
    case2 = df1[df1[:,group] .== cases[2], yvar]
    test = MannWhitneyUTest(case1,case2)
end

function survivalrate_algorythm(var; step = 0.05, xaxis = nothing)
    isnothing(xaxis) && (xaxis = range(extrema(var)..., step = step))
    survival = 1 .- ecdf(var).(xaxis)
    return (Xaxis = collect(xaxis), fy = survival)
end

function cumulative_algorythm(var; step = 0.05, xaxis = nothing)
    isnothing(xaxis) && (xaxis = range(extrema(var)..., step = step))
    cum = ecdf(var).(xaxis)
    return (Xaxis = collect(xaxis), fy = cum)
end

function hazardrate_algorythm(var; step = 0.05, xaxis = nothing)
    isnothing(xaxis) && (xaxis = range(extrema(var)..., step = step))
    survival = 1 .- ecdf(var).(xaxis)
    hazard = -pushfirst!(diff(survival),0)./survival
    return (Xaxis = collect(xaxis), fy = hazard)
end
