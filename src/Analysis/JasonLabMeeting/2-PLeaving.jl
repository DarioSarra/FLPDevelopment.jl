############################# Age
Age_LevingAn, Age_LevingAn_df = function_analysis(Age_s,:LogDuration, cumulative_algorythm; 
    grouping = :Age, calc = :bootstrapping)
xprop = ("Poke Time(seconds)",(log10.([0.1,1,10,100,1000]),
    ["0.1","1","10","100","1000"]))
yprop = ("Probablity of leaving")
plot(Age_LevingAn, xaxis = xprop, yaxis = yprop, legend = false,
    yticks = 0:0.2:1,
    left_margin = -8mm,
    bottom_margin = -8mm,
    right_margin = 1.5mm,
    tickfont = font(10, "Bookman Light"))
savefig(joinpath(temp_dir,"Age_pleaving.pdf"))
############################# Sex
##
Sex_LevingAn, Sex_LevingAn_df = function_analysis(Age_s,:LogDuration, cumulative_algorythm; 
    grouping = :Sex, calc = :bootstrapping)
xprop = ("Poke Time(seconds)",(log10.([0.1,1,10,100,1000]),
    ["0.1","1","10","100","1000"]))
yprop = ("Probablity of leaving")
plot(Sex_LevingAn, legend = false,
    yticks = 0:0.2:1,
    left_margin = -8mm,
    bottom_margin = -8mm,
    right_margin = 1.5mm,
    color = [:violetred3 :sandybrown],
    tickfont = font(10, "Bookman Light"))
savefig(joinpath(temp_dir,"Sex_pleaving.pdf"))
############################# Virus
##
Cas_LevingAn, Cas_LevingAn_df = function_analysis(Cas_s,:LogDuration, cumulative_algorythm; 
    grouping = :Virus, calc = :bootstrapping)
xprop = ("Poke Time(seconds)",(log10.([0.1,1,10,100,1000]),
    ["0.1","1","10","100","1000"]))
yprop = ("Probablity of leaving")
plot(Cas_LevingAn, xaxis = xprop, yaxis = yprop, legend = false,
    yticks = 0:0.2:1,
    left_margin = -8mm,
    bottom_margin = -8mm,
    right_margin = 1.5mm,
    tickfont = font(10, "Bookman Light"))
savefig(joinpath(temp_dir,"Cas_pleaving.pdf"))