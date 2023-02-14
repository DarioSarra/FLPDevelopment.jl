include("0-PlotSetting.jl")
############################# Age
##
Age_Rewards = Difference(Age_s, :Age, :Num_Rewards; ylims = (0,1.15),
    ylabel = "Number of rewards per trial", color = [1,2],
    yticks = 0:0.2:1.5)
    Age_Rewards.plt
Age_Rewards.test
Age_Rewards.groupdf
savefig(joinpath(temp_dir,"age_nrew.pdf"))
############################# Virus
##
Cas_Rewards = Difference(Cas_s, :Virus, :Num_Rewards;
    font_size = 6,
    ylims = (0,1.20),
    ylabel = "Number of rewards per trial", color = [1,2],
    yticks = 0:0.2:1.5)
    Cas_Rewards.plt
Cas_Rewards.test
Cas_Rewards.groupdf
savefig(joinpath(temp_dir,"Cas_nrew.pdf"))