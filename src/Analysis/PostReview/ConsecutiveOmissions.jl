##
Age_AfterLast = Difference(Age_s, :Age, :AfterLast,
    ylabel = "Consecutive failures before leaving", ylims = (0,5), size = (350,600),
    guidefont = font(18, "Bookman Light"),xtickfont = font(18, "Bookman Light"))
    Age_AfterLast.plt
Age_AfterLast.test
Age_AfterLast.groupdf
Age_AfterLast.groupdf[!,:Measure] .= "ConsecutiveFailures"
rename!(Age_AfterLast.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_ConsecutiveFailures.png"))
##
##
Cas_AfterLast = Difference(Cas_s, :Virus, :AfterLast,
    ylabel = "Consecutive failures before leaving", ylims = (0,5), size = (350,600),
    guidefont = font(18, "Bookman Light"),xtickfont = font(18, "Bookman Light"))
    Cas_AfterLast.plt
Cas_AfterLast.test
Cas_AfterLast.groupdf
Cas_AfterLast.groupdf[!,:Measure] .= "ConsecutiveFailures"
rename!(Cas_AfterLast.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_ConsecutiveFailures.png"))
