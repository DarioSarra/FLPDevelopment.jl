include(joinpath(pwd(),"src","Analysis","Filters.jl"))
##
Age_Rewards = Difference(Age_s, :Age, :Num_Rewards; ylims = (0,1.15),
    ylabel = "number of rewards per trial")
    Age_Rewards.plt
Age_Rewards.test
Age_Rewards.groupdf
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Review", "Rewards.pdf"))
##
Cas_Rewards = Difference(Cas_s, :Virus, :Num_Rewards; ylims = (0,1.23),
    ylabel = "Number of rewards per trial")
    Cas_Rewards.plt
Cas_Rewards.test
Cas_Rewards.groupdf
