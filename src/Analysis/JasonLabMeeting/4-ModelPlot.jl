################################ Age
dtype = [String, Union{Missing, String}, String, String, Float64, String, String]
Age_bdf = CSV.read(joinpath(temp_dir, "Age_Leaving_1000Bootstrap.csv"), DataFrame; types = dtype)
Age_bdf.err = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in Age_bdf.err]
@df Age_bdf scatter(:coef, :Variable, xerror = :err, 
    legend = false, markercolor = :grey,
    markersize = 7,
    left_margin = -13mm,
    bottom_margin = -9mm,
    xtickfont = 6,
    ytickfont = 8,
    guidefont = font(10, "Bookman Light"),
    formatter = :plain,
    xticks = (-1:0.2:0.7,string.(-1:0.2:0.7)),
    xlabel = "Estimated coefficients \n (1000 samp)")
    vline!([0,0], color = :red, linestyle = :dash)
savefig(joinpath(temp_dir, "Age_Bootstrap.pdf"))
################################ Sex
## simple
dtype = [String, Union{Missing, String}, String, String, Float64, String, String]
Sex_bdf = CSV.read(joinpath(temp_dir, "Sex_simple_1000Bootstrap.csv"), DataFrame; types = dtype)
Sex_bdf.err = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in Sex_bdf.err]
@df Sex_bdf scatter(:coef, :Variable, xerror = :err, 
    legend = false, markercolor = :grey,
    markersize = 7,
    left_margin = -13mm,
    bottom_margin = -9mm,
    xtickfont = 6,
    ytickfont = 8,
    guidefont = font(10, "Bookman Light"),
    formatter = :plain,
    xticks = (-1:0.2:0.8,string.(-1:0.2:0.8)),
    xlabel = "Estimated coefficients \n (1000 samp)")
    vline!([0,0], color = :red, linestyle = :dash)
savefig(joinpath(temp_dir, "Sex_simple_Bootstrap.pdf"))
## full

@df Sex_complete_bdf scatter(:coef, 1:1:length(:Variable), xerror = :err, 
    legend = false, markercolor = :grey,
    markersize = 5,
    size = (400,500),
    left_margin = 0mm,
    bottom_margin = -9mm,
    xtickfont = 4,
    ytickfont = 3,
    yticks = (collect(1:1:length(:Variable)),:Variable),
    guidefont = font(10, "Bookman Light"),
    formatter = :plain,
    xticks = (-1.4:0.2:1,string.(-1.4:0.2:1)),
    xlabel = "Estimated coefficients \n (1000 samp)")
    vline!([0,0], color = :red, linestyle = :dash)
savefig(joinpath(temp_dir, "Sex_complete_Bootstrap.pdf"))

################################ Virus
##
dtype = [String, Union{Missing, String}, String, String, Float64, String, String]
Cas_bdf = CSV.read(joinpath(temp_dir, "Cas_Leaving_1000Bootstrap.csv"), DataFrame; types = dtype)
Cas_bdf.err = [tuple(parse.(Float64,split(x[2:end-1], ", "))...) for x in Cas_bdf.err]
@df Cas_bdf scatter(:coef, :Variable, xerror = :err, 
    legend = false, markercolor = :grey,
    markersize = 7,
    left_margin = -13mm,
    bottom_margin = -9mm,
    xtickfont = 5,
    ytickfont = 8,
    guidefont = font(10, "Bookman Light"),
    formatter = :plain,
    xticks = (-1.2:0.2:1,string.(-1.2:0.2:1)),
    xlabel = "Estimated coefficients \n (1000 samp)")
    vline!([0,0], color = :red, linestyle = :dash)
savefig(joinpath(temp_dir, "Cas_Bootstrap.pdf"))