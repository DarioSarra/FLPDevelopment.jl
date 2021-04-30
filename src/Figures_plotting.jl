function add_bar!(plt, ref, span)
    plot!(plt,[1,1],[ref,ref+span], linecolor = :black)
    plot!(plt,[2,2],[ref,ref+span], linecolor = :black)
    plot!(plt,[1,2],[ref+span/2,ref+span/2], linecolor = :black)
end

function add_pvalue!(plt, ref, span, p)
    if isa(p, String)
        message = p
    else
        message = p < 0.01 ? "p < 0.01" : "p < 0.05"
    end
    annotate!(plt,[(1.5,ref + span,
        Plots.text(message,
        6, :center))])
    return plt
end

function add_effect!(plt, ref, span, e)
    message = "Effect size = $e"
    annotate!(plt,[(1.5,ref - span,
        Plots.text(message,
        6, :center))])
    return plt
end

function add_info!(plt, df, p, e; ylabel = "Median ...")
    plt_lim = maximum(df.Central .+ last.(df.ERR))
    plt_span = plt_lim/10
    ref, span = plt_lim + plt_span, plt_span/2
    add_bar!(plt, ref, span)
    if p >= 0.05
        pval = "N.S."
    else
        round(p, digits = 2) == 0 ? pval = round(p, digits = 3) : pval = round(p, digits = 2)
        add_effect!(plt, ref+span/2, span/2, e)
    end
    pmessage = "Mann-Whitney U test, p = $(pval)"
    add_pvalue!(plt, ref+span/2, span/2, pmessage)
    yaxis!(ylabel = ylabel, ylims = (0,plt_lim + 2plt_span))
    xlabel!("Group")
    return plt
end

function Difference(df, group, var; ind_summary = mean, ylabel = "Median ...")
    res_plt, res_group, res_individual = median_ci_scatter(df, group, var; ind_summary = ind_summary)
    res_test = FLPDevelopment.MWU_test(res_individual, group, var)
    res_effect = res_test.U/(res_test.nx * res_test.ny)
    if res_effect < 0.5
        res_effect = ((res_test.nx * res_test.ny) - res_test.U) / (res_test.nx * res_test.ny)
    end
    e = round(res_effect; digits =2)
    add_info!(res_plt, res_group, pvalue(res_test), e; ylabel = ylabel)
    return (plt = res_plt,
        test = res_test,
        effect = res_effect,
        groupdf = res_group,
        individual_df = res_individual)
end
