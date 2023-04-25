include(joinpath(pwd(),"src","Analysis","Filters.jl"))
##
Age_Npokes = Difference(Age_s, :Age, :Num_pokes, ylabel = "Number of pokes per trial", ylims = (0,5.5))
Age_Npokes.plt
Age_Npokes.test
Age_Npokes.groupdf
Age_Npokes.groupdf[!,:Measure] .= "Number of pokes"
rename!(Age_Npokes.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Age_NumPokes_test.pdf"))
##
Age_AfterLast = Difference(Age_s, :Age, :AfterLast, ylabel = "Consecutive pokes after last reward", ylims = (0,5.5))
Age_AfterLast.plt
Age_AfterLast.test
Age_AfterLast.groupdf
Age_AfterLast.groupdf[!,:Measure] .= "Number of pokes"
rename!(Age_AfterLast.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Age_AfterLast_test.pdf"))
##
Cas_Npokes = Difference(Cas_s, :Virus, :Num_pokes, ylabel = "Number of pokes per trial", ylims = (0,5.5))
Cas_Npokes.plt
Cas_Npokes.test
Cas_Npokes.groupdf
Cas_Npokes.groupdf[!,:Measure] .= "Number of pokes"
# rename!(Cas_Npokes.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Cas_NumPokes_test.pdf"))
##
Cas_AfterLast = Difference(Cas_s, :Virus, :AfterLast, ylabel = "Consecutive pokes after last reward", ylims = (0,5.5))
Cas_AfterLast.plt
Cas_AfterLast.test
Cas_AfterLast.groupdf
Cas_AfterLast.groupdf[!,:Measure] .= "Number of pokes"
# rename!(Cas_AfterLast.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","April2023","Cas_AfterLast_test.pdf"))
