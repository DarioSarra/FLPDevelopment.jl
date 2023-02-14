count_correct(v) = sum(v[2:end]) / (length(v)-1)
count_incorrect(v) = 1 - count_correct(v)
################################ Age
Age_corr1 = combine(groupby(Age_p,[:MouseID,:Streak,:Age]), :Correct => first => :CorrectStart)
Age_corr2 = combine(groupby(Age_corr1,[:MouseID,:Age]), :CorrectStart => sum, nrow => :Trials)

Age_Inc = Difference(Age_corr1, :Age, :CorrectStart; 
    ind_summary = count_correct, font_size = 6, 
    ylabel = "Fraction of correct leaving", ylims = (0,1.25),
    yticks = 0:0.2:1,
    size = (360,500),
    left_margin = -5.5mm,
    bottom_margin = -6mm,
    widen = false)
Age_Inc.plt
Age_Inc.groupdf
Age_Inc.test
savefig(joinpath(replace(path,basename(path)=>""),"Development_Figures","Fig1","Incorrect_test.pdf"))
################################ Virus
##
Cas_corr1 = combine(groupby(Cas_p,[:MouseID,:Streak,:Virus]), :Correct => first => :CorrectStart)
Cas_corr2 = combine(groupby(Cas_corr1,[:MouseID,:Virus]), :CorrectStart => sum, nrow => :Trials)

Cas_Inc = Difference(Cas_corr1, :Virus, :CorrectStart; 
    ind_summary = count_incorrect, font_size = 6, 
    ylabel = "Fraction of correct leaving", ylims = (0,1.25),
    yticks = 0:0.2:1,
    size = (360,500),
    left_margin = -5.5mm,
    bottom_margin = -6mm,
    widen = false)
Cas_Inc.plt
Cas_Inc.groupdf
Cas_Inc.test