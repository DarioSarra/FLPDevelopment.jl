include(joinpath(splitpath(@__DIR__)[1:end-1]...,"Filters.jl"))
##
println(propertynames(Age_p))
Age_gdf = combine(groupby(Age_p, [:MouseID,:Age]), [:Reward,:PokeOut] => 
((r,t) -> cumsum(r)./t) => :Rewards_over_time)
Age_RewardRate = Difference(Age_gdf, :Age, :Rewards_over_time; ylims = (0,0.05),
    ylabel = "Reward rate per second")
    Age_RewardRate.plt
Age_RewardRate.test
Age_RewardRate.groupdf
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023", "Age_RewardRate.png"))
##
Cas_gdf = combine(groupby(Cas_p, [:MouseID,:Virus]), [:Reward,:PokeOut] => 
((r,t) -> cumsum(r)./t) => :Rewards_over_time)
Cas_RewardRate = Difference(Cas_gdf, :Virus, :Rewards_over_time; ylims = (0,0.05),
    ylabel = "Reward rate per second")
    Cas_RewardRate.plt
Cas_RewardRate.test
Cas_RewardRate.groupdf
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023", "Cas_RewardRate.png"))
