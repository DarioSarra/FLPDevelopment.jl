################################ Age
Age_Npokes = Difference(Age_s, :Age, :Num_pokes, font_size = 6,
    ylabel = "Number of pokes per trial", ylims = (0,5.5),
    # size = (360,500),
    color = [1,2],
    left_margin = -7mm,
    bottom_margin = -6mm,
    widen = false)
Age_Npokes.plt
Age_Npokes.test
Age_Npokes.groupdf
Age_Npokes.groupdf[!,:Measure] .= "Number of pokes"
rename!(Age_Npokes.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(temp_dir,"Age_NumPokes.pdf"))
################################ Virus
##
Cas_Npokes = Difference(Cas_s, :Virus, :Num_pokes, font_size = 6,
    ylabel = "Number of pokes per trial", ylims = (0,6.5),
    size = (360,500),
    color = [1,2],
    left_margin = -5.5mm,
    bottom_margin = -6mm,
    widen = false)
Cas_Npokes.plt
Cas_Npokes.test
Cas_Npokes.groupdf
Cas_Npokes.groupdf[!,:Measure] .= "Number of pokes"
rename!(Cas_Npokes.groupdf,[:Central => :Median, :ERR => :CI])
savefig(joinpath(temp_dir,"Cas_NumPokes.pdf"))