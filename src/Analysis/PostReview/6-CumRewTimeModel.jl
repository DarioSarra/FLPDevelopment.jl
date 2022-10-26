include(joinpath(splitpath(@__DIR__)[1:end-1]...,"Filters.jl"))
##############################
Age_CR_verb = @formula(Leave ~ 1 + Streak_zscore*Age + LogOut_zscore*Age + CumReward*Age +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Age_CR = fit(MixedModel,Age_CR_verb, Age_p, Bernoulli(); contrasts)
AIC_test(Age_CR,Age_Full)
aic(Age_CR)
aic(Age_Full)
##
# dtype = [String, Union{Missing, String}, String, String, Float64, String, String]
# Age_CR_bdf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_CumRew_Bootstrap.csv"), DataFrame; types = dtype)
# Age_CR_bdf = CSV.read(joinpath("/home/beatriz/Documents/Datasets/Development_figures", "Age_CumRew_1000Bootstrap.csv"), DataFrame; types = dtype)
# Age_CR_bdf.err = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in Age_CR_bdf_2.err]
Age_CR_bdf = bootstrapdf(Age_p,Age_CR; n = 1000)
transform!(Age_CR_bdf,:names => ByRow(x -> replace(x,r"centered: \d{1}" => "",
    ")" => "", "("=> "",
    "Streak_zscore" => "Trial",
    "LogOut_zscore" => "PokeTime",
    "&" => "&\n")) => :Variable)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_CumRew_1000Bootstrap.csv"),Age_CR_bdf)
CSV.write(joinpath("/home/beatriz/Documents/Datasets/Development_figures", "Age_CumRew_1000Bootstrap.csv"),Age_CR_bdf)
@df Age_CR_bdf scatter(:coef, :Variable, xerror = :err, legend = false, markercolor = :grey)
    vline!([0,0], color = :red, linestyle = :dash, xlabel = "Estimated coefficients (1000 samp)")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Age_CumRew_Bootstrap.png"))
##
Cas_CR_verb = @formula(Leave ~ 1 + Streak_zscore*Virus + LogOut_zscore*Virus + CumReward*Virus +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Cas_CR = fit(MixedModel,Cas_CR_verb, Cas_p, Bernoulli(); contrasts)
AIC_test(Cas_CR,Cas_Full)
aic(Cas_CR)
aic(Cas_Full)
##
# Cas_CR_bdf = CSV.read(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_CumRew_Bootstrap.csv"), DataFrame; types = dtype)
Cas_CR_bdf = CSV.read(joinpath("/home/beatriz/Documents/Datasets/Development_figures", "Cas_CumRew_1000Bootstrap.csv"), DataFrame; types = dtype)
Cas_CR_bdf.err = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in Cas_CR_bdf_2.err]
Cas_CR_bdf = bootstrapdf(Cas_p,Cas_CR; n = 1000)
transform!(Cas_CR_bdf,:names => ByRow(x -> replace(x,r"centered: \d{1}" => "",
    ")" => "", "("=> "",
    "Streak_zscore" => "Trial",
    "LogOut_zscore" => "PokeTime",
    "&" => "&\n")) => :Variable)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_CumRew_1000Bootstrap.csv"),Cas_CR_bdf)
CSV.write(joinpath("/home/beatriz/Documents/Datasets/Development_figures", "Cas_CumRew_1000Bootstrap.csv"),Cas_CR_bdf)
@df Cas_CR_bdf scatter(:coef, :Variable, xerror = :err, legend = false, markercolor = :grey)
    vline!([0,0], color = :red, linestyle = :dash, xlabel = "Estimated coefficients (1000 samp)")
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Cas_CumRew_Bootstrap.png"))
##
