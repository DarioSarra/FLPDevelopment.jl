include(joinpath(pwd(),"src","Analysis","Filters.jl"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","Age_PokeData.csv"), select(Age_p,[:MouseID,:Age,:Sex,:Side,:Streak,:Poke,:Correct,:In, :Out]))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","Virus_PokeData.csv"), select(Cas_p,[:MouseID,:Virus,:Sex,:Side,:Streak,:Poke,:Correct,:In, :Out]))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","Age_TrialData.csv"), select(Age_s,[:MouseID,:Age,:Sex,:Side,:Streak,:Num_pokes,:LogDuration]))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","Virus_TrialData.csv"), select(Cas_s,[:MouseID,:Virus,:Sex,:Side,:Streak,:Num_pokes,:LogDuration]))

##
Age_Npokes = Difference(Age_s, :Age, :Num_pokes, ylabel = "Number of pokes per trial", ylims = (0,5.5))
Age_Npokes.plt
Age_Npokes.test
Age_Npokes.groupdf
Age_Npokes.groupdf[!,:Measure] .= "Number of pokes"
rename!(Age_Npokes.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Age_NumPokes_test.pdf"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","Fig2","F_Age_Npokes_test.csv"), select(Age_Npokes.individual_df, Not([:xpos])))
##
Age_AfterLast = Difference(Age_s, :Age, :AfterLast, ylabel = "Consecutive pokes after last reward", ylims = (0,5.5))
Age_AfterLast.plt
Age_AfterLast.test
Age_AfterLast.groupdf
Age_AfterLast.groupdf[!,:Measure] .= "Number of pokes"
rename!(Age_AfterLast.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Age_AfterLast_test.pdf"))c
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","Fig2","F_Age_AfterLast_test.csv"), select(Age_AfterLast.individual_df, Not([:xpos])))
##
Cas_Npokes = Difference(Cas_s, :Virus, :Num_pokes, ylabel = "Number of pokes per trial", ylims = (0,5.5))
Cas_Npokes.plt
Cas_Npokes.test
Cas_Npokes.groupdf
Cas_Npokes.groupdf[!,:Measure] .= "Number of pokes"
# rename!(Cas_Npokes.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Cas_NumPokes_test.pdf"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","Fig2","F_Cas_Npokes_test.csv"), select(Cas_Npokes.individual_df, Not([:xpos])))

##
Cas_AfterLast = Difference(Cas_s, :Virus, :AfterLast, ylabel = "Consecutive pokes after last reward", ylims = (0,5.5))
Cas_AfterLast.plt
Cas_AfterLast.test
Cas_AfterLast.groupdf
Cas_AfterLast.groupdf[!,:Measure] .= "Number of pokes"
# rename!(Cas_AfterLast.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Cas_AfterLast_test.pdf"))
CSV.write(joinpath(replace(path,basename(path)=>""),"Development_Figures","Submission2023","Fig3","F_Cas_AfterLast_test.csv"), select(Cas_AfterLast.individual_df, Not([:xpos])))
