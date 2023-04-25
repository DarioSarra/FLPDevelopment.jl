include(joinpath(pwd(),"src","Analysis","Filters.jl"))
##
Age_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Age_Basic = fit(MixedModel,Age_Basic_verb, Age_p, Bernoulli();contrasts)
Age_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Age + LogOut_zscore * Age +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Age_Full = fit(MixedModel,Age_Full_verb, Age_p, Bernoulli();contrasts)
MixedModels.likelihoodratiotest(Age_Basic,Age_Full)
Age_noTrial_verb = @formula(Leave ~ 1 + LogOut_zscore * Age +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Age_noTrial =  fit(MixedModel,Age_noTrial_verb, Age_p, Bernoulli();contrasts)
MixedModels.likelihoodratiotest(Age_noTrial,Age_Full)

##
Age_BootDf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","1000AgeBootstrap.csv"),DataFrame)
# Age_BootDf = bootstrapdf(Age_p, Age_Full; n = 1000)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","1000AgeBootstrap.csv"),Age_BootDf)
Age_btdf = Age_BootDf[[1,4,2,3,5,6],:]
open_html_table(Age_btdf)
Age_btdf.err = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in Age_btdf.err]
# Age_btdf.variable
yprop = ("",font(8, "Bookman Light"),(collect(1:nrow(Age_btdf)),
    ["Intercept","PokeTime","Trial","Juveniles",
    "Trial &\nJuveniles",
    "PokeTime &\nJuveniles"]))
    xprop = ("Coefficient estimate", font(8, "Bookman Light"), (-1.33,1.3))
    @df Age_btdf scatter(:coef ,1:nrow(Age_btdf), xerror = :err,legend = false,
    xaxis = xprop, yaxis = yprop, markercolor = :gray75, leftmargin = -10mm)
    vline!([0], linecolor = :red, legend = false, linestyle = :dash)
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review","Age_LevingRateBootstrap.pdf"))
##
Cas_Basic_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Cas_Basic = fit(MixedModel,Cas_Basic_verb, Cas_p, Bernoulli();contrasts)
Cas_Full_verb = @formula(Leave ~ 1 + Streak_zscore * Virus + LogOut_zscore * Virus +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Cas_Full = fit(MixedModel,Cas_Full_verb, Cas_p, Bernoulli();contrasts)
MixedModels.likelihoodratiotest(Cas_Basic,Cas_Full)
##
Cas_noTrial_verb = @formula(Leave ~ 1 + LogOut_zscore * Virus +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID));
Cas_noTrial =  fit(MixedModel,Cas_noTrial_verb, Cas_p, Bernoulli();contrasts)
MixedModels.likelihoodratiotest(Cas_noTrial,Cas_Full)
