Sex_LevingAn, Sex_LevingAn_df = function_analysis(Age_s,:LogDuration, cumulative_algorythm;
    grouping = :Sex, calc = :bootstrapping, color = [:violetred3 :sandybrown])
    xprop = ("Poke Time(seconds)", xyfont,(log10.([0.1,1,10,100,1000]),["0.1","1","10","100","1000"]))
    yprop = ("Probablity of leaving", xyfont)
    plot!(Sex_LevingAn, xaxis = xprop, yaxis = yprop, legend = (0.8, 0.25))
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","SexLevingRate.pdf"))
rename!(Sex_LevingAn_df, [:Central => :Median, :ERR => :CI, :LogDuration => :PokeTime])
sort!(Sex_LevingAn_df,:Sex)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission","SFig3", "C","SexLevingRate.csv"), select(Sex_LevingAn_df, Not([:up,:low]))
##
Age_Sex_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + CumReward + Sex*Age +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Age_Sex = fit(MixedModel,Age_Sex_verb, Age_p, Bernoulli(); contrasts)
AIC_test(Age_Sex,Age_CR)
##
Age_Sex_bdf = bootstrapdf(Age_p,Age_Sex; n = 100)

transform!(Age_Sex_bdf,:names => ByRow(x -> replace(x,r"centered: \d{1}" => "",
    ")" => "", "("=> "",
    "Streak_zscore" => "Trial",
    "LogOut_zscore" => "PokeTime",
    "&" => "&\n")) => :Variable)
@df Age_Sex_bdf scatter(:coef, :Variable, xerror = :err, legend = false, markercolor = :grey)
    vline!([0,0], color = :red, linestyle = :dash, xlabel = "Estimated coefficients (100 samp)")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Sex_CumRew_Bootstrap.pdf"))

##
Age_Sex0_verb = @formula(Leave ~ 1 + Streak_zscore*Age + LogOut_zscore*Age + CumReward*Age + Sex*Age +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Age_Sex0 = fit(MixedModel,Age_Sex0_verb, Age_p, Bernoulli(); contrasts)
AIC_test(Age_Sex0,Age_CR)
