
Age_p.Occupancy_zscore

Age_occupancy_base_form = @formula(Occupancy_zscore ~ 1 + LogOut_zscore + (1+LogOut_zscore|MouseID));
Age_occupancy_full_form = @formula(Occupancy_zscore ~ 1 + LogOut_zscore * Age + (1+LogOut_zscore|MouseID));
Age_Occupancy_base = fit(MixedModel,Age_occupancy_base_form, Age_p)
Age_occupancy_full = fit(MixedModel,Age_occupancy_full_form, Age_p)
MixedModels.likelihoodratiotest(Age_Occupancy_base ,Age_occupancy_full)

OccupancyAge_BootDf = bootstrapdf(occupancy_df, Occupancy_m2; n = 1000)
# OccupancyAge_BootDf.names = ["Intercept", "PokeTime","Group:Juveniles","PokeTime&Group:Juveniles"]
# OccupancyAge_BootDf.Y = 1:nrow(OccupancyAge_BootDf)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","1000AgeOccupancyBootstrap.csv"),OccupancyAge_BootDf)
##
Cas_p.Occupancy_zscore

Cas_occupancy_base_form = @formula(Occupancy_zscore ~ 1 + LogOut_zscore + (1+LogOut_zscore|MouseID));
Cas_occupancy_full_form = @formula(Occupancy_zscore ~ 1 + LogOut_zscore * Virus + (1+LogOut_zscore|MouseID));
Cas_Occupancy_base = fit(MixedModel,Cas_occupancy_base_form, Cas_p)
Cas_occupancy_full = fit(MixedModel,Cas_occupancy_full_form, Cas_p)
MixedModels.likelihoodratiotest(Cas_Occupancy_base ,Cas_occupancy_full)

OccupancyCas_BootDf = bootstrapdf(occupancy_df, Occupancy_m2; n = 1000)
# OccupancyCas_BootDf.names = ["Intercept", "PokeTime","Group:Juveniles","PokeTime&Group:Juveniles"]
# OccupancyCas_BootDf.Y = 1:nrow(OccupancyCas_BootDf)
# CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","1000AgeOccupancyBootstrap.csv"),OccupancyCas_BootDf)