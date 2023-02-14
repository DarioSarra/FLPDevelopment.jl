############################# Age
## Calculate Model
Age_base_formula = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + CumReward +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Age_base = fit(MixedModel,Age_base_formula, Age_p, Bernoulli(); contrasts)

Age_Leaveing_formula = @formula(Leave ~ 1 + Streak_zscore*Age + LogOut_zscore*Age + CumReward*Age +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Age_Leaving = fit(MixedModel,Age_Leaveing_formula, Age_p, Bernoulli(); contrasts)

Age_test = MixedModels.likelihoodratiotest(Age_base,Age_Leaving)
## Calculate bootstrap
# Age_bdf = bootstrapdf(Age_p,Age_Leaving; n = 1000)
# transform!(Age_bdf,:names => ByRow(x -> replace(x,r"centered: \d{1}" => "",
#     ")" => "", "("=> "",
#     "Streak_zscore" => "Trial",
#     "LogOut_zscore" => "PokeTime",
#     "&" => "&\n")) => :Variable)

# Age_bdf = Age_bdf[[1,4,5,2,3,7,8,6],:]
# CSV.write(joinpath(temp_dir, "Age_Leaving_1000Bootstrap.csv"),Age_bdf)
## Load Bootstrap
dtype = [String, Union{Missing, String}, String, String, Float64, String, String]
Age_bdf = CSV.read(joinpath(temp_dir, "Age_Leaving_1000Bootstrap.csv")) => DataFrame
Age_bdf.err = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in Age_bdf.err]
################################ Sex
## Calculate model 1
Age_Sex_verb = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + CumReward + Sex*Age +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Age_Sex = fit(MixedModel,Age_Sex_verb, Age_p, Bernoulli(); contrasts)
Sex_test = MixedModels.likelihoodratiotest(Age_base,Age_Sex)
##Calculate bootstrap
Age_Sex_bdf = bootstrapdf(Age_p,Age_Sex; n = 1000)
transform!(Age_Sex_bdf,:names => ByRow(x -> replace(x,r"centered: \d{1}" => "",
    ")" => "", "("=> "",
    "Streak_zscore" => "Trial",
    "LogOut_zscore" => "PokeTime",
    "&" => "&\n")) => :Variable)
Age_Sex_bdf = Age_Sex_bdf[[1,3,4,2,6,5,7],:]
CSV.write(joinpath(temp_dir, "Sex_simple_1000Bootstrap.csv"),Age_Sex_bdf)
## Calculate model 2
Age_Sex_verb2 = @formula(Leave ~ 1 + Streak_zscore*Age*Sex + LogOut_zscore*Age*Sex + CumReward*Age*Sex +
(1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Age_Sex2 = fit(MixedModel,Age_Sex_verb2, Age_p, Bernoulli(); contrasts)
Sex_test = MixedModels.likelihoodratiotest(Age_Leaving,Age_Sex2)
##Calculate bootstrap
Age_Sex_bdf2 = bootstrapdf(Age_p,Age_Sex2; n = 1000)
transform!(Age_Sex_bdf2,:names => ByRow(x -> replace(x,r"centered: \d{1}" => "",
    ")" => "", "("=> "",
    "Streak_zscore" => "Trial",
    "LogOut_zscore" => "PokeTime",
    "&" => "&\n")) => :Variable)
Sex_complete_bdf = Age_Sex_bdf2[[1,5,6,2,3,10,12,7,4,11,13,8,9,15,16,14],:]
CSV.write(joinpath(temp_dir, "Sex_complete_1000Bootstrap.csv"),Age_Sex_bdf2)
################################ Virus
## Calculate model
Cas_base_formula = @formula(Leave ~ 1 + Streak_zscore + LogOut_zscore + CumReward +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Cas_base = fit(MixedModel,Cas_base_formula, Cas_p, Bernoulli(); contrasts)
Cas_Leaveing_formula = @formula(Leave ~ 1 + Streak_zscore*Virus + LogOut_zscore*Virus + CumReward*Virus +
    (1|MouseID)+(Streak_zscore|MouseID)+(LogOut_zscore|MouseID)+(CumReward|MouseID));
Cas_Leaving = fit(MixedModel,Cas_Leaveing_formula, Cas_p, Bernoulli(); contrasts)
Cas_test = MixedModels.likelihoodratiotest(Cas_base,Cas_Leaving)
## Calculate bootstrap
# Cas_bdf = bootstrapdf(Cas_p,Cas_Leaving; n = 1000)
# transform!(Cas_bdf,:names => ByRow(x -> replace(x,r"centered: \d{1}" => "",
#     ")" => "", "("=> "",
#     "Streak_zscore" => "Trial",
#     "LogOut_zscore" => "PokeTime",
#     "&" => "&\n")) => :Variable)
# Cas_bdf = Cas_bdf[[1,4,5,2,3,7,8,6],:]
# CSV.write(joinpath(temp_dir, "Cas_Leaving_1000Bootstrap.csv"),Cas_bdf)
## load bootstrap
dtype = [String, Union{Missing, String}, String, String, Float64, String, String]
Cas_bdf = CSV.read(joinpath(temp_dir, "Cas_Leaving_1000Bootstrap.csv"),DataFrame)
Cas_bdf.err = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in Cas_bdf.err]
