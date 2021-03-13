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
function individual_ecdf(df,xvar; err = :MouseID, axis = nothing)
    if isnothing(axis)
        axis = ecdf(df[:,xvar]).sorted_values
    end
    gd1 = groupby(df, err)
    df1 = combine(gd1) do dd
        CUM = ecdf(dd[:,xvar])
        (Xaxis = collect(axis), Vals = CUM(axis))
    end
end

function individual_cdf(df,xvar; err = :MouseID, axis = nothing)

    gd1 = groupby(df, err)
    df1 = combine(gd1) do dd
        df2 = DataFrame(Xaxis = Int[], Vals = Float64[])
        for (k,i) in countmap(dd[:,xvar])
            push!(df2, [k,i/nrow(dd)])
        end
        df2
    end
    sort!(df1,:Xaxis)
    gdp = groupby(df1,err)
    transform!(gdp,:Vals => cumsum => :Vals)
    return sort(df1, :Xaxis)
end

function group_cdf(df,xvar; err = :MouseID, group = nothing, estimated = false)
    axis = minimum(df[:,xvar]):maximum(df[:,xvar])
    if estimated
        cumulative = individual_ecdf
    else
        cumulative = individual_cdf
    end
    if isnothing(group)
        df = cumulative(df,xvar; err = err, axis = axis)
        gd = groupby(df,:Xaxis)
        return combine(gd, :Vals => mean => :Mean, :Vals => sem => :Sem)
    else
        gd1 = groupby(df, group)
        res = combine(gd1) do dd
            m_res = cumulative(dd,xvar; err = err, axis = axis)
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

group_distribution(df, var; args...) = length(unique(df[:,var])) > 30 ? group_kde(df, var; args ...) : group_frequency(df, var; args ...)

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
    #evaluate if it's going over multiple grouping columns
    multiple = typeof(xvar) <: AbstractVector && sizeof(xvar) > 1
    if multiple
        gdc = groupby(df,vcat(xvar,err))
    else
        gdc = groupby(df,[xvar,err])
    end
    df1 = combine(yvar => summary => yvar,gdc)
    sort!(df1,xvar)
    if multiple
        label_df = unique(df1[:,xvar])
    else
        label_df = unique(df1[:,[xvar]])
    end
    label_df[!,:xpos] = 1:nrow(label_df)
    leftjoin(df1,label_df; on = xvar)
    firstval = union(df1[:,xvar])[1]
    df1[!,:xpos] = [v == firstval ? 1 : 2  for v in df1[:,xvar]]
    return df1
end

function group_summary(df1,xvar,yvar; normality = true)
    if normality
        central = mean
    else
        central = median
    end
    df2 = combine(groupby(df1,xvar)) do dd
        # using bootstrap to calculate the central point and the 95% CI
        b = bootstrap(central,dd[:,yvar], BasicSampling(100000))
        #bootstrap confint returns a tuple containing a tuple with the statistic
        #estimates and the CI lower and uppper bounds
        m, lower, upper = confint(b,BasicConfInt(0.95))[1]
        ci1 = m - lower
        ci2 = upper - m
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
        legend = false,
        markersize = 10,
        markerstroke = 1)
    if showmice
        @df df1 scatter!(:xpos,cols(yvar),
        markersize = 5,
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

function overtrial_plot(df, grouping, var)
    #calculate mean of var for each mouse over a trial's bin
    dfSum = summary_xy(df,:BinnedStreak,var; summary = mean, err = :MouseID, group = grouping)
    #plot mean over mice plus sem
    overtrialplot = @df dfSum plot(:BinnedStreak,:Mean, ribbon = :Sem, group = cols(grouping), linecolor = :auto,
        legend = :topright, ylabel = String(var), xlabel = "Trial")
    return overtrialplot
end

function frequency_plot(df, grouping, var)
    df1 = group_frequency(df,var; group = grouping)
    filter!(r -> !isnan(r.Sem), df1)
    @df df1 plot(:Xaxis, :Mean, ribbon = :Sem, group = cols(grouping), linecolor = :auto)
end

function cdf_plot(df, grouping, var; estimated = false)
    df1 = group_cdf(df,var; group = grouping, estimated = estimated)
    filter!(r -> !isnan(r.Sem), df1)
    @df df1 plot(:Xaxis, :Mean, ribbon = :Sem, group = cols(grouping), linecolor = :auto,
     legend=:bottomright, xlabel = String(var), ylabel = "Cumulative distribution")
end

function mean_sem_scatter(df, grouping, var)
    df1 = individual_summary(df, grouping, var)
    df2 = combine(groupby(df1,grouping)) do dd
        m = mean(dd[:,var])
        s = sem(dd[:,var])
        (Central = m, ERR = (s,s))
    end
    if typeof(grouping) <: AbstractVector && sizeof(grouping) > 1
        df2[!,:xaxis] = [join(x,"_") for x in eachrow(df2[:, grouping])]
        xaxis = :xaxis
    else
        xaxis = grouping
    end
    plt = @df df2 scatter(1:nrow(df2),:Central,
        yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),cols(xaxis)),
        legend = false)
    return plt, df2
end

function mode_ci_scatter(df, grouping, var)
    df1 = individual_summary(df, grouping, var)
    df2 = group_summary(df1, grouping, var; normality = false)
    if typeof(grouping) <: AbstractVector && sizeof(grouping) > 1
        df2[!,:xaxis] = [join(x,"_") for x in eachrow(df2[:, grouping])]
        xaxis = :xaxis
    else
        xaxis = :grouping
    end
    plt = @df df2 scatter(1:nrow(df2),:Central,
        yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),cols(xaxis)),
        legend = false)
    return plt, df2
end

function mean_ci_scatter(df, grouping, var)
    df1 = individual_summary(df, grouping, var)
    df2 = group_summary(df1, grouping, var; normality = true)
    @df df2 scatter(1:nrow(df2),:Central,
        yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),cols(grouping)),
        legend = false)
end

function incorrect_fraction_scatter(df, grouping, var)
    df1 = individual_summary(df, grouping, var; summary = incorrect_fraction)
    df2 = combine(groupby(df1,grouping)) do dd
        m = mean(dd[:,var])
        s = sem(dd[:,var])
        (Central = m, ERR = (s,s))
    end
    if typeof(grouping) <: AbstractVector && sizeof(grouping) > 1
        df2[!,:xaxis] = [join(x,"_") for x in eachrow(df2[:, grouping])]
        xaxis = :xaxis
    else
        xaxis = grouping
    end
    plt = @df df2 scatter(1:nrow(df2),:Central,
        yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),cols(xaxis)),
        legend = false)
    return plt, df2
end

function var_plots(df,var; grouping = grouping)
    # bins along streaks to plot over trials
    df[!,:BinnedStreak] = bin_axis(df.Streak; unit_step = 5)
    overtrialplot = overtrial_plot(df,var; grouping = grouping)

    #prepare violin plot, each group density and dotplot has to be calculate separetely
    vals = union(df[:,grouping])
    vplot = @df  df[df[:,grouping] .== vals[1],:] violin(string.(cols(grouping)),cols(var),
         label = vals[1], side = :left, legend = :top)
     @df df[df[:,grouping] .== vals[2],:] violin!(string.(cols(grouping)),cols(var),
        label = vals[2], side = :right)
    @df df[df[:,grouping] .== vals[1],:] dotplot!(string.(cols(grouping)),cols(var),
        label = "", side = :left, legend = :top, markersize = 4, markercolor = :grey, markeralpha = 0.7)
    @df df[df[:,grouping] .== vals[2],:] dotplot!(string.(cols(grouping)),cols(var),
        label = "", side = :right, legend = :top, markersize = 4, markercolor = :grey,
        markeralpha = 0.7)

    # scatter plot of mean of the mice's mean and sem
    df1 = individual_summary(df, grouping, var)
    df2 = group_summary(df1, grouping, var; normality = true)
    statisticplot = @df df2 scatter(1:nrow(df2),:Central,
        yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),cols(grouping)),
        legend = false)

    return overtrialplot, vplot, statisticplot
end

function vars_fig(df_s; grouping = nothing, summary_opt = :MEAN)
    if isnothing(grouping)
        if in(:Age,propertynames(df_s))
            grouping = :Age
        elseif in(:Virus,propertynames(df_s))
            grouping = :Virus
         end
    end
    # return variable mean over trials, distribution, and statistics
    ALt,ALd,ALs = var_plots(df_s, :AfterLast; grouping = grouping)
    ILt,ILd,ILs = var_plots(df_s, :IncorrectLeave; grouping = grouping)
    plot(ALt,ALd,ALs,
        ILt,ILd,ILs,
        layout = grid(2,3),
        size=(1200,1200),
        thickness_scaling = 1)
end

# Function to add a title over multiple figures
function maintitle!(plt,what::String)
    y = ones(3)
    title = Plots.scatter(y, marker=0, markeralpha=0, grid = false,
    axis=false, leg=false,ticks = nothing, size = (200,1200),
    annotations=(2, y[2], Plots.text(what)))
    plot(title, plt, layout = grid(2,1, heights = [0.05,0.95]))
end

function median_duration(df, group, variable)
    if group == :Age
        firstval = "Adults"
        vals = ["Adults", "Juveniles"]
    elseif group == :Virus
        firstval = "tdTomato"
        vals = ["tdTomato", "Caspase"]
    end
    gd = groupby(df, [:MouseID, group])
    df2 = combine(gd, variable => median => variable)
    df2.xpos = [x == firstval ? 1 : 2 for x in df2[:, group]]
    plt = @df df2 scatter(:xpos, cols(variable),
        xlims = (0.5, 2 + 0.5),
        # xticks = (:xpos,cols(group)),
        xticks = ([1,2], vals),
        legend = false)
    return plt, df2
end

function duration_analysis(df, group, variable)
    # med_scat, med_df = FLPDevelopment.median_duration(df, group, variable)
    median_analysis = DoubleAnalysis(df,group, variable;summary_opt = :MEDIAN, showmice = true)
    xlabel!(median_analysis.nonparametric_plot, String(group))
    ylabel!(median_analysis.nonparametric_plot, String(variable))
    k = group_kde(df,variable, group = group)
    dist = @df k plot(:Xaxis,:Mean, ribbon = :Sem, group = cols(group), linecolor = :auto,
    ylabel = "PDF", xlabel = String(variable))
    # return dist, med_scat, med_df
    return dist, median_analysis.nonparametric_plot
end

function all_duration_analysis(poke, streak, group)
    trial_d, trial_s = duration_analysis(streak,group, :Trial_duration)
    pokep = filter(r -> r.PreInterpoke > 0, poke)
    inter_d, inter_s = duration_analysis(pokep,group, :PreInterpoke)
    poke_d, poke_s = duration_analysis(pokep,group, :PokeDur)
    travel_d, travel_s = duration_analysis(streak,group, :Travel_to)

    p1 = plot(trial_d, travel_d,
        trial_s, travel_s,
        layout = grid(2,2),
        size=(874,620),
        thickness_scaling = 1)

    p2 = plot(inter_d,poke_d,
        inter_s, poke_s,
        layout = grid(2,2),
        size=(874,620),
        thickness_scaling = 1)
    return p1,p2
end

function bin_duration(df,variable, group; modality = :QUANTILE, qs = 0.1:0.1:0.9, length = 10, unit_step = nothing, fontx = 11)
    mice_bin = calculate_bin_duration(df,variable, group; modality = modality, qs = qs, length = length, unit_step = unit_step)
    group_bin = combine(groupby(mice_bin,:Bin), :Count=> mean, :Count => sem)
    return plot_bin_duration(mice_bin, group; variable = String(variable), fontx = fontx), group_bin
end

function calculate_bin_duration(df,variable, group; modality = :QUANTILE, qs = 0.1:0.1:0.9, length = 10, unit_step = nothing)
    df1 = copy(df[:, [:MouseID, group, variable]])
    future_col = Symbol("Binned_" * String(variable))
    if modality == :QUANTILE
        treshvec = quantile_bin(df1[:,variable]; qs = qs)
    elseif modality == :STEP
        treshvec = range_bin(df1[:,variable]; length = length, unit_step = unit_step)
    end
    binvec = [findfirst(x .<= treshvec) for x in df1[:,variable]]
    valvec = [isnothing(x) ? treshvec[end] : treshvec[x] for x in binvec]
    df1[!, future_col] = valvec
    gd = groupby(df1,[:MouseID, group])
    df3 = combine(gd) do dd
        df2 = DataFrame(Bin = Float64[], Count = Int[])
        # for (k,v) in countmap(dd[:, future_col])
        for t in treshvec
            push!(df2,[t,sum(dd[:,future_col] .== t)])
        end
        df2
    end
    sort!(df3,[:Bin,:MouseID])
    # rename!(df3, :Bin => Symbol("Binned_" * String(variable)))
end

function quantile_bin(vec; qs = 0.1:0.1:0.9)
    dist = fit(Gamma,vec)
    treshvec = quantile.(dist,qs)
    push!(treshvec, treshvec[end] + treshvec[end] - treshvec[end-1])
    return treshvec
end

function range_bin(v; length = 10, unit_step = nothing)
    if isnothing(unit_step)
        treshvec = collect(range(extrema(v)...;length = length))
        push!(treshvec, treshvec[end] + treshvec[end] - treshvec[end-1])
        return treshvec
    else
        treshvec = collect(range(extrema(v)...;step = unit_step))
        push!(treshvec, treshvec[end] + treshvec[end] - treshvec[end-1])
        return treshvec
    end
end

function plot_bin_duration(df, group; variable = "undef", fontx = 11)
    plots = []
    for m in union(df.MouseID)
        dd = filter(r -> r.MouseID == m, df)
        ntrials = sum(dd.Count)
        if group == :Virus
            c = dd[1,:Virus] == "Caspase" ? :red : :black
        elseif group == :Age
            c = dd[1,:Age] == Juveniles ? :red : :black
        else
            c = :auto
        end
        plt = @df dd bar(:Bin,:Count,
            xticks = round.(union(:Bin), digits = 2),
            xrotation = 45,
            xtickfontsize = fontx,
            xlabel = String(variable),
            yticks = 0:5:100,
            grid = true,
            color = c,
            ylabel = m,
            label = "n trials = " * string(ntrials),
            legend = :top)
        push!(plots, plt)
    end
    return plots
end
