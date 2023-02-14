include(joinpath(splitpath(@__DIR__)[1:end-1]...,"Filters.jl"))
gr(aspect_ratio = :none,
    left_margin = -9mm,
    size = (500,500),
    bottom_margin = -5mm,
    thickness_scaling = 1.5,
    guidefont = font(14, "Bookman Light"),
    tickfont = font(12, "Bookman Light"))
    temp_dir = "/Users/dariosarra/Documents/Lab/Oxford/Lerch/Presentations/Lab_meeting/20230210"