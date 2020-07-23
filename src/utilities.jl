function dvplot(df,xvar,yvar,stat = :parametric)
    gdc = groupby(df,[xvar,:MouseID])
    df1 = combine(yvar => mean => yvar,gdc)
    sort!(df1,xvar)
    firstval = union(df1[:,xvar])[1]
    df1[!,:xpos] = [v == firstval ? 1 : 2  for v in df1[:,xvar]]
    df2 = combine(groupby(df1,xvar)) do dd
        m = mean(dd[:,yvar])
        ci1 = m - confint(OneSampleTTest(dd[:,yvar]))[1]
        ci2 = confint(OneSampleTTest(dd[:,yvar]))[2] - m
        (Mean = mean(dd[:,yvar]), ERR = (ci1,ci2))
    end
    @df df2 scatter(1:nrow(df2),:Mean, yerror = :ERR,
        xlims = (0.5, nrow(df2) + 0.5),
        xticks = (1:nrow(df2),cols(xvar)),
        legend = false)
    @df df1 scatter!(:xpos,cols(yvar), markersize = 3, alpha = 0.5, color = :grey)
end
