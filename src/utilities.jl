function summarize(dd::AbstractDataFrame,Xvar::Symbol,Yvar::Symbol; Err = :MouseID , mode = :sem)
    ErrGroups = vcat(Xvar,Err)
    XaxisGroups = vcat(Xvar)
    pre_err = combine(groupby(dd, ErrGroups)) do df
        (Mean = mean(df[:,Yvar]),)
    end
    if mode == :sem
        with_err = combine(groupby(pre_err,XaxisGroups)) do df
            (Mean = mean(df.Mean), SEM = sem(df.Mean))
        end
        rename!(with_err, Xvar=>:Xaxis)
        sort!(with_err,:Xaxis)
        filter!(r -> !isnan(r.SEM), with_err)
    elseif mode == :conf
        with_err = combine(groupby(pre_err,XaxisGroups)) do df
            ci = confint(OneSampleTTest(df.Mean))
            m = mean(df.Mean)
            (Mean = m, ERRlow = m - ci[1], ERRup = ci[2] - m)
        end
        with_err[!,:ERR] = [(low,up) for (low,up) in zip(with_err.ERRlow,with_err.ERRup)]
        rename!(with_err, Xvar=>:Xaxis)
        sort!(with_err,:Xaxis)
    end
    return with_err
end
