# function NormalDifference(df, group, var; ind_summary = mean, ylabel = "Median ...", xyfont = font(18, "Bookman Light"))
#     res_plt, res_group, res_individual = mean_sem_scatter(df, group, var; ind_summary = ind_summary)
#     res_test = test_difference(res_individual, group, var; normality = true)
#     pooledSD = sqrt(((res_group.Central[2])^2 + (res_group.Central[2])^2)/2)
#     res_effect = (res_group.Central[2] - res_group.Central[2]) / pooledSD
#     # if res_effect < 0.5
#     #     res_effect = ((res_test.nx * res_test.ny) - res_test.U) / (res_test.nx * res_test.ny)
#     # end
#     e = round(res_effect; digits =2)
#     add_info!(res_plt, res_group, pvalue(res_test), e; normality = true)
#     xprop = ("Group", xyfont)
#     yprop = (ylabel, xyfont)
#     plot!(xaxis = xprop, yaxis = yprop)
#     return (plt = res_plt,
#         test = res_test,
#         effect = res_effect,
#         groupdf = res_group,
#         individual_df = res_individual)
# end

"""
    `Difference(df, group, var; ind_summary = mean, ylabel = "Median ...", xyfont = font(18, "Bookman Light"), ylims = nothing)`

    A - plot a scatter with median ± CI of the variable var per each value of group
    B - perform MannWhitneyUTest and annotate the results
    C - Calculate the effect size if p < 0.05
    D - return a named tuple with the plot, MannWhitneyUTest, effect size
        a Dataframe with the group median and CI of var,
        and a Dataframe with the individual mouse mean value of var
"""
function Difference(df, group, var; ind_summary = mean, font_size = 8, kwargs...)
    res_plt, res_group, res_individual = median_ci_scatter(df, group, var; ind_summary = ind_summary, kwargs...)
    res_test = MWU_test(res_individual, group, var)
    res_effect = res_test.U/(res_test.nx * res_test.ny)
    if res_effect < 0.5
        res_effect = ((res_test.nx * res_test.ny) - res_test.U) / (res_test.nx * res_test.ny)
    end
    e = round(res_effect; digits = 2)
    add_info!(res_plt, res_group, pvalue(res_test), e; font_size = font_size)
    # xprop = ("Group", xyfont)
    # yprop = (ylabel, xyfont)
    # plot!(xaxis = xprop, yaxis = yprop)
    return (plt = res_plt,
        test = res_test,
        effect = res_effect,
        groupdf = res_group,
        individual_df = res_individual)
end

"""
    `median_ci_scatter(df, grouping, var)`
    A - calculate the mean of the variable var for each animal using individual summary
    B - calculate the median and ci for each group using group_summary
    C - Plots the value as a scatter plot ± CI
"""

function median_ci_scatter(df, grouping, var; ind_summary = mean, kwargs...)
    df1 = individual_summary(df, grouping, var; summary = ind_summary)
    df2 = group_summary(df1, grouping, var; normality = false)
    if typeof(grouping) <: AbstractVector && sizeof(grouping) > 1
        df2[!,:xaxis] = [join(x,"_") for x in eachrow(df2[:, grouping])]
        xaxis = :xaxis
    else
        xaxis = grouping
    end
    try
        cases = levels(df[:,grouping])
        df2[!,:color] = [val == cases[1] ? OGCol1 : OGCol2 for val in df2[:,xaxis]]
    catch
        df2[!,:color] .= :gray75
    end
    plt = @df df2 scatter(1:nrow(df2),:Central,
        yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),string.(cols(xaxis))),
        color = :color,
        legend = false; kwargs...)
    return plt, df2, df1
end

function add_info!(plt, df, p, e; normality = false, font_size = 8)
    ylims = Plots.ylims(plt)
    plt_lim = isnothing(ylims) ? maximum(df.Central .+ last.(df.ERR)) : last(ylims) - (last(ylims) - first(ylims))/5
    plt_bot = isnothing(ylims) ? minimum(df.Central .+ first.(df.ERR)) : first(ylims)
    plt_span = (plt_lim-plt_bot)/10
    ref, span = plt_lim + plt_span, plt_span/2
    add_bar!(plt, ref, span)
    if p >= 0.05
        pval = "N.S."
    else
        round(p, digits = 2) == 0 ? pval = round(p, digits = 3) : pval = round(p, digits = 2)
        add_effect!(plt, ref+span/2, span/2, e; font_size = font_size)
    end
    if normality
        pmessage = "T-test, p = $(pval)"
    else
        pmessage = "Mann-Whitney U test, p = $(pval)"
    end
    add_pvalue!(plt, ref+span/2, span, pmessage; font_size = font_size)
    isnothing(ylims) && (ylims = (0,plt_lim + 2plt_span))
    yaxis!(ylims = ylims)
    return plt
end

function add_bar!(plt, ref, span)
    plot!(plt,[1,1],[ref,ref+span], linecolor = :black)
    plot!(plt,[2,2],[ref,ref+span], linecolor = :black)
    plot!(plt,[1,2],[ref+span/2,ref+span/2], linecolor = :black)
end

function add_pvalue!(plt, ref, span, p; font_size = 8)
    if isa(p, String)
        message = p
    else
        message = p < 0.01 ? "p < 0.01" : "p < 0.05"
    end
    annotate!(plt,[(1.5,ref + span,
        Plots.text(message,
        font_size, :center))])
    return plt
end

function add_effect!(plt, ref, span, e; font_size = 8)
    message = "Effect size = $e"
    annotate!(plt,[(1.5,ref - span,
        Plots.text(message,
        font_size, :center))])
    return plt
end

"""
    `function_analysis(df,variable, f; grouping = nothing, step =0.05, calc = :basic,
            color = [:auto], linestyle = [:auto])`

    Apply the function f over the variable var per each value of grouping and
    plots the result over the variable var. functions algorythms are defined in Analysis_fun.jl
"""

function function_analysis(df,var, f; grouping = nothing, step =0.05, calc = :basic,
        color = [:auto], linestyle = [:auto])
    subgroups = isnothing(grouping) ? [:MouseID] : vcat(:MouseID,grouping)
    xaxis = range(extrema(df[:, var])..., step = step)
    dd1 = combine(groupby(df,subgroups), var => (t-> f(t,xaxis = xaxis)) => AsTable)
    rename!(dd1, Dict(:Xaxis => var))
    sort!(dd1,[:MouseID,var])
    if calc == :bootstrapping
        dd2 = combine(groupby(dd1,grouping)) do dd3
            group_summary(dd3,var,:fy; normality = false)
        end
        dd2[!,:low] = [x[1] for x in dd2.ERR]
        dd2[!,:up] = [x[2] for x in dd2.ERR]
    elseif calc == :quantiles
        dd2 = combine(groupby(dd1,[grouping,var]), :fy =>(t-> (Central = mean(t),
        low= abs(mean(t) - quantile(t,0.25),
        up = abs(quantile(t,0.975)-mean(t))),
        # ERR = (abs(mean(t) - quantile(t,0.25)) + abs(quantile(t,0.975)-mean(t)))/2,
        SEM = sem(t))) => AsTable)
    elseif calc == :basic
        dd2 = combine(groupby(dd1,[grouping,var]), :fy =>(t-> (Central = mean(t),up = sem(t), low = sem(t))) => AsTable)
    end
    sort!(dd2,var)

    plt = @df dd2 plot(cols(var),:Central, ribbon = (:low, :up), group = cols(grouping), linecolor = :auto, color = color, linestyle = linestyle)
    return plt, dd2
end

function mediansurvival_analysis(streakdf,variable, grouping; plt = plot())
    dd1 = combine(groupby(streakdf,[:MouseID,grouping]), variable => median => variable)
    # dd2 = combine(groupby(dd1,grouping), variable => (t-> (Mean = mean(t),Sem = sem(t))) => AsTable)
    dd2 = group_summary(dd1,grouping,variable; normality = false)
    # dd2[!,:low] = [x[1] for x in dd2.ERR]
    # dd2[!,:up] = [x[2] for x in dd2.ERR]
    isa(dd2[:, grouping], CategoricalArray) && (dd2[!,grouping] = string.(dd2[!,grouping]))
    println(isa(dd2[:, grouping], CategoricalArray))
    @df dd2 scatter!(plt,cols(grouping), :Central, yerror = :ERR,
        # xlims = (-0.25,2.25), xlabel = "Group",
        ylabel = "Median survival time", label = "")
    return plt
end
##
function plot_xy(df,x,y; group = "none", bin = true, digits = 1, calc = :mean, filtNaN = true, kwargs...)
    res = summarize_xy(df,x,y; group = group, bin = bin, digits = digits, calc = calc)
    bin ? (res_x = Symbol("Binned"*string(x))) : (res_x = x)
    y_m = Symbol(string(y) * "_mean")
    y_s = Symbol(string(y) * "_sem")
    filtNaN && filter!(y_s => v -> !isnan(v),res)
    if group == "none"
        return @df res plot(cols(res_x),cols(y_m), ribbon = cols(y_s); kwargs...)
    else
        return @df res plot(cols(res_x),cols(y_m), ribbon = cols(y_s), group = cols(group); kwargs...)
    end
end

function summarize_xy(df,x,y; group = "none", bin = true, digits = 1, calc = :mean)
    if bin
        xcol = Symbol("Binned"*string(x))
        df[:,xcol] = round.(df[:,x],digits = digits)
    else
        xcol = x
    end
    standardgroup = [xcol]
    group == "none" ? (grouping = standardgroup) : (grouping = vcat(standardgroup,group))
    df1 = combine(groupby(df,vcat(grouping,:MouseID)), y => mean => y)
    if calc == :mean
        df2 = combine(groupby(df1,grouping), y .=> [mean, sem])
    elseif calc == :bootstrapping
        df2 = combine(groupby(df, group)) do dd
            group_summary(dd,xcol,y; normality = false)#
        end
        df2[!,:low] = [x[1] for x in df2.ERR]
        df2[!,:up] = [x[2] for x in df2.ERR]
    end
    sort!(df2,grouping)
end
"""
    `pdf_heatmap(df,group,var)`
    Calculate PDF of 'var' per mouse split by 'group' and returns a heatmap over 100 bins of 'var' 
"""
function pdf_heatmap(df,group,var)
    g = vcat(group,:MouseID)
    df0 = groupby(df,g)
    df1 = combine(df0, var => kde => :KDE)
    r = range(extrema(df[:,var])..., length = 100)
    r = range(-1,3, length = 100)
    df2 = transform(df1,:KDE => ByRow(k-> map(c->pdf(k,c), r)) =>:PDF)
    df3 = sort(df2[:,Not(:KDE)],group)
    pos_dic = countmap(df3[:,group])
    samples = map(x->pos_dic[x], levels(df3[:,group]))
    stepping(x1,x2) = x1*2 +x2/2
    pos_ticks = accumulate(stepping,samples, init = 0)
    # return df3, pos_dic, pos_ticks
    PDFmat = reduce(hcat,df3.PDF)'
    heatmap(PDFmat,yticks = (pos_ticks, string.(levels(df3[:,group]))),
        xticks = (collect(range(0,100,length = 5)), ["0","1","10","100","1000"]),
        xlabel = "Leaving time density (s)",
        color= :tol_ylorbr
    )
    hline!(cumsum(collect(values(pos_dic))).+0.5, color = :black, label ="")
end