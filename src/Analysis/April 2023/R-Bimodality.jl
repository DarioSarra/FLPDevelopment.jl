include(joinpath(splitpath(@__DIR__)[1:end-1]...,"Filters.jl"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Age_s.csv"),Age_s[:,[:Age,:MouseID,:LogDuration,:Strategy]])
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Cas_s.csv"),Cas_s[:,[:Virus,:MouseID,:LogDuration,:Strategy]])
## Cluster LogDuration data in 2 kmeans
kcluster_df!(Age_s,2,:LogDuration)
categorise_kclusters!(Age_s, :LogDuration, :Assignment, ["Short", "Long"]; new_name = :Strategy)
## calculate proportion of strategies per mouse
Age_k = combine(groupby(Age_s,[:Age,:MouseID])) do dd
    diz = countmap(dd.Strategy)
    (Short = diz["Short"]/nrow(dd), Long = diz["Long"]/nrow(dd), Trials = nrow(dd))
end
#calculate difference, but per se we don't expect the difference to be equal to 0
transform!(Age_k, [:Long,:Short] => ByRow((l,s)-> l-s) => :Diff)
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Age_k.csv"),Age_k)

kcluster_df!(Cas_s,2,:LogDuration)
categorise_kclusters!(Cas_s, :LogDuration, :Assignment, ["Short", "Long"]; new_name = :Strategy)
## calculate proportion of strategies per mouse
Cas_k = combine(groupby(Cas_s,[:Virus,:MouseID])) do dd
    diz = countmap(dd.Strategy)
    (Short = diz["Short"]/nrow(dd), Long = diz["Long"]/nrow(dd), Trials = nrow(dd))
end
#calculate difference, but per se we don't expect the difference to be equal to 0
transform!(Cas_k, [:Long,:Short] => ByRow((l,s)-> l-s) => :Diff)
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Cas_k.csv"),Cas_k)

