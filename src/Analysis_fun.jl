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

"""
    `bootstrapdf(df, mdl; grouping = nothing, n = 100)`
    bootstrap a model results and saves it in a dataframes to store as csv
"""
function bootstrapdf(df, mdl; grouping = nothing, n = 100)
    rng = MersenneTwister(1234321)
    samp1 = parametricbootstrap(rng,n,mdl)
    sampdf = DataFrame(samp1.allpars)
    bootdf = combine(groupby(sampdf,[:type, :group, :names]), :value => shortestcovint => :interval)
    filter!(r -> ismissing(r.group), bootdf)
    bootdf.coef = coef(mdl)
    isnothing(grouping) && (grouping = check_group(df))
    cases = levels(df[:,grouping])
    try (bootdf.variable = ["Intercept", "Trial", "$grouping: $(cases[2])", "Poke-time",
        "Trial & $grouping: $(cases[2])",
        "Poke-time & $grouping: $(cases[2])"])
    catch
        println("grouping renaming didn't work")
    end
    transform!(bootdf, [:coef, :interval] => ByRow((c,e) -> (c -e[1], e[2]-c)) => :err)
    return bootdf
end

"""
    `peaksdf(df,var;group = nothing)`
    performs a KDE of the vector 'var, finds the peaks
    and saves it in a dataframes to store as csv
"""

function peaksdf(df,var;group = nothing, flat = true)
    isnothing(group) ? (df0 = df) : (df0 = groupby(df,group))
    df1 = combine(df0, var => kde => :KDE)
    transform!(df1, :KDE => ByRow(y -> findmaxima(y.density)) => [:Pos,:Dens])
    transform!(df1,[:KDE, :Pos] => ByRow((y,p) -> y.x[p]) => :Val)
    if flat
        df2 = flatten(df1,[:Pos,:Dens,:Val])[:,Not(:KDE)]
    else
        df2 = transform(df1,:Pos => ByRow(length) => :Count)[:,Not(:KDE)]
    end
    return df2
end

"""
    `outliersdf(df,var;group = nothing)`
    1-performs a KDE on the vector 'var' by mouse and group, 
    2-finds the each KDE peaks and saves it in a dataframes
    3-calclulate the overall KDE of 'var' not split by group
    4-for each peak returns the p value that it belongs to the 'var' distribution
    5-in case of multiple peaks returns the joint probability that they blong to the 'var' distribution
    6-verify if any animal has a joint probability lower than 5%
"""

function outliersdf(df,var;group = nothing)
    g = isnothing(group) ? [:MouseID] : vcat(group,:MouseID)
    micepeaks = FLPDevelopment.peaksdf(df,var; group = g, flat = false)
    k_density = kde(df[:,var])
    transform!(micepeaks, :Val => ByRow(x -> 
        map(c-> 
        pdf(k_density,c)#=/(sum(pdf(k_density,range(-1,3,length=100))))=#,
        x)) 
        => :p)
    joint_p(v) = length(v) == 1 ? v[1] : sum(v) - prod(v)
    # transform!(micepeaks, :p => ByRow(x->sum(x) - *(vcat(1.00,x)...)) => :p_sum)
    transform!(micepeaks, :p => ByRow(joint_p) => :p_union)
    transform!(micepeaks, :p_union => ByRow(x -> x<0.05) => :Outliers)
    return micepeaks
end
